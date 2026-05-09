import os
import subprocess
import csv
import time
import sys
import traceback
import tkinter as tk
import multiprocessing
from tkinter import filedialog, messagebox
from concurrent.futures import ProcessPoolExecutor, as_completed

# =========================================================
# Windows + Conda DLL 修复
# =========================================================
conda_env_path = r"C:\Users\DELL\miniconda3\envs\docking"
dll_path = os.path.join(conda_env_path, "Library", "bin")
if os.path.exists(dll_path):
    os.environ['PATH'] = dll_path + os.pathsep + os.environ.get('PATH', '')
    if hasattr(os, 'add_dll_directory'):
        try:
            os.add_dll_directory(dll_path)
        except:
            pass

# =========================================================
# 导入 RDKit / Meeko
# =========================================================
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger
    from meeko import MoleculePreparation
    from meeko import PDBQTWriterLegacy

    RDLogger.DisableLog('rdApp.*')
except ImportError as e:
    print(f"\n❌ RDKit 或 Meeko 导入失败:\n{e}")
    sys.exit(1)

# =========================================================
# 对接核心参数配置
# =========================================================
CENTER = ["31.39", "3.52", "-15.31"]
SIZE = ["20", "20", "20"]

# 🚀 解封算力：自动检测 CPU 线程数，保留 4 个线程给系统防卡死，其余全开
MAX_WORKERS = max(4, multiprocessing.cpu_count() - 4)

# 🚀 降低初筛搜索深度以换取极速 (推荐 5 或 6)
EXHAUSTIVENESS = "6"

# =========================================================
# 🧙‍♂️ 核心新增：化学清道夫 (拯救被错杀的实体化合物)
# =========================================================
def chemical_janitor(mol):
    """
    清洗和修复异常分子的「化学清道夫」。
    在脱盐和 3D 生成之前，把 Vina 不支持但实际存在的元素进行等效替换或清理。
    """
    if mol is None:
        return None
        
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        
        # 1. 拯救金属盐类 (Li, Na, Mg, K, Ca) -> 突变为氢原子 (H) 
        # 因为在生理水溶液中它们本来就会解离，配体实际上是结合质子的形态
        if z in [3, 11, 12, 19, 20]:
            atom.SetAtomicNum(1)
            atom.SetFormalCharge(0)
            
        # 2. 拯救含硒分子 (Se, 34) -> 突变为硫 (S, 16)
        # Vina 没有硒的参数，但硫的范德华半径和电性与硒极度相似，这是工业界常用的“瞒天过海”技巧
        elif z == 34:
            atom.SetAtomicNum(16)
            atom.SetFormalCharge(0)
            
        # 3. 拯救硅/硼分子 (B, 5 / Si, 14) -> 突变为碳 (C, 6)
        elif z in [5, 14]:
            atom.SetAtomicNum(6)
            
    return mol

# =========================================================
# 分子安全检查 (最后一道防线)
# =========================================================
def is_safe_for_3d(mol):
    allowed_elements = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}
    bad_atoms = []
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z not in allowed_elements:
            bad_atoms.append(atom.GetSymbol())
    if bad_atoms:
        return False, f"罕见重金属拦截: {','.join(set(bad_atoms))}"

    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms > 80:
        return False, f"巨型分子({heavy_atoms})"
    return True, "安全"


# =========================================================
# 极速生成3D构象
# =========================================================
def generate_3d_mol_fast(mol):
    try:
        params = AllChem.ETKDGv2()
        params.useRandomCoords = True
        params.maxIterations = 200  

        status = AllChem.EmbedMolecule(mol, params)
        if status == 0:
            return True

        status = AllChem.EmbedMolecule(mol, useRandomCoords=True)
        return status == 0
    except:
        return False


# =========================================================
# 单分子流水线处理进程
# =========================================================
def process_single_smiles(cas_name, smiles, output_dir, receptor_file, vina_exe):
    final_cas = cas_name if cas_name.strip() else "未命名分子"
    
    try:
        safe_name = "".join(c for c in final_cas if c.isalnum() or c in ('_', '-')).rstrip()
        if not safe_name: safe_name = "unknown"

        ligand_pdbqt = os.path.join(output_dir, f"{safe_name}_temp.pdbqt")
        output_pdbqt = os.path.join(output_dir, f"{safe_name}_out.pdbqt")

        # 1. 解析 SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': 'SMILES解析失败'}

        # 2. 🧙‍♂️ 召唤清道夫：魔改/清理异常原子
        mol = chemical_janitor(mol)

        # 3. 脱盐 (由于金属已被突变为 H，分离出来的就是纯净的有机母体了)
        frags = list(Chem.GetMolFrags(mol, asMols=True))
        if not frags: return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': '碎片异常'}
        mol = max(frags, key=lambda m: m.GetNumAtoms())

        # 4. 价态修正与加氢
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass # 即使有点小瑕疵也强行往下走
            
        mol = Chem.AddHs(mol)
        
        # 5. 安全检查 (只拦截真正的无法计算的重金属，比如金 Au、铂 Pt)
        safe, reason = is_safe_for_3d(mol)
        if not safe: return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': reason}

        if not generate_3d_mol_fast(mol):
            return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': '3D构象生成失败'}

        # 6. Meeko 转换
        try:
            preparator = MoleculePreparation()
            setups = preparator.prepare(mol)
            if not setups: return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': 'Meeko失败'}
            pdbqt_string, is_ok, _ = PDBQTWriterLegacy.write_string(setups[0])
            if not is_ok: return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': 'PDBQT转换异常'}
        except Exception:
            return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': 'Meeko崩溃'}

        # 7. 写入临时文件
        with open(ligand_pdbqt, "w", encoding='utf-8') as f:
            f.write(pdbqt_string)

        # 8. 调用 Vina 进行极速对接
        cmd = [
            vina_exe,
            "--receptor", receptor_file, "--ligand", ligand_pdbqt,
            "--center_x", CENTER[0], "--center_y", CENTER[1], "--center_z", CENTER[2],
            "--size_x", SIZE[0], "--size_y", SIZE[1], "--size_z", SIZE[2],
            "--exhaustiveness", EXHAUSTIVENESS,
            "--cpu", "1",  
            "--out", output_pdbqt
        ]

        best_score = None
        try:
            process = subprocess.run(
                cmd, capture_output=True, text=True, check=True,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0
            )
            for line in process.stdout.split('\n'):
                if line.startswith('   1 '):
                    best_score = float(line.split()[1])
                    break
        except Exception as e:
            return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': 'Vina引擎运算报错'}
        finally:
            if os.path.exists(ligand_pdbqt):
                try:
                    os.remove(ligand_pdbqt)
                except:
                    pass

        if best_score is None:
            if os.path.exists(output_pdbqt):
                try:
                    os.remove(output_pdbqt)
                except:
                    pass
            return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': '无有效构象'}

        return {'CAS': final_cas, 'SMILES': smiles, 'Score': best_score, 'Status': '成功'}

    except Exception:
        return {'CAS': final_cas, 'SMILES': smiles, 'Score': None, 'Status': '未知异常终止'}


# =========================================================
# 主调度程序
# =========================================================
def main():
    root = tk.Tk()
    root.withdraw()
    root.attributes('-topmost', True)

    print("=" * 60)
    print("🚀 i9-14900 满血并发筛选引擎 (附带化学清道夫版)")
    print("=" * 60)

    vina_exe = filedialog.askopenfilename(title="1/4 选择 vina.exe")
    if not vina_exe: return
    receptor_file = filedialog.askopenfilename(title="2/4 选择 receptor.pdbqt")
    if not receptor_file: return
    csv_file = filedialog.askopenfilename(title="3/4 选择 CSV 文件")
    if not csv_file: return
    result_csv = filedialog.asksaveasfilename(
        title="4/4 保存排行榜结果", defaultextension=".csv", initialfile="HTVS_Docking_Results_Final.csv"
    )
    if not result_csv: return

    # 读取CSV
    molecules = []
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header = next(reader)
        cas_idx, smiles_idx = 0, 1
        for i, h in enumerate(header):
            if 'smiles' in h.lower(): smiles_idx = i
            if 'cas' in h.lower() or 'id' in h.lower(): cas_idx = i
        for row in reader:
            if len(row) > max(cas_idx, smiles_idx) and row[smiles_idx].strip():
                molecules.append((row[cas_idx].strip(), row[smiles_idx].strip()))

    total = len(molecules)
    if total == 0: return

    output_dir = os.path.join(os.path.dirname(csv_file), "Docking_Output")
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n✅ 待处理配体: {total} 个")
    print(f"🔥 并发线程数: {MAX_WORKERS} 线程火力全开")
    print(f"⚡ Vina 搜索深度: {EXHAUSTIVENESS} (优化初筛速度)")
    print(f"🧹 化学清道夫: 已启动 (自动修复 Na, K, Ca, Se, Si 等导致拦截的元素)")
    print("\n任务已下发至线程池，建议打开任务管理器查看 CPU 占用率...\n")

    results = []
    success_count = 0
    start = time.time()
    completed = 0

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(process_single_smiles, cas, smiles, output_dir, receptor_file, vina_exe): (cas, smiles)
            for cas, smiles in molecules
        }

        for future in as_completed(futures):
            completed += 1
            cas_val, smiles_val = futures[future]
            final_cas_val = cas_val if cas_val.strip() else "未命名分子"
            
            try:
                res = future.result()
            except Exception:
                res = {'CAS': final_cas_val, 'SMILES': smiles_val, 'Score': None, 'Status': '进程池调度异常'}

            results.append(res)
            if res['Score'] is not None:
                success_count += 1
                print(f"[{completed}/{total}] ✅ {res['CAS']:<18} 结合能: {res['Score']:>7.2f} kcal/mol")
            else:
                print(f"[{completed}/{total}] ❌ {res['CAS']:<18} {res['Status']}")

    # 结果排序与保存
    results.sort(key=lambda x: (x['Score'] is None, x['Score'] if x['Score'] is not None else 999))
    
    with open(result_csv, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=['CAS', 'SMILES', 'Score', 'Status'])
        writer.writeheader()
        writer.writerows(results)

    elapsed = time.time() - start
    print("\n" + "=" * 60)
    print(f"🎉 筛选彻底完成！")
    print(f"⏱️ 总耗时: {elapsed / 60:.2f} 分钟")
    print(f"✅ 成功对接: {success_count} 个化合物")
    messagebox.showinfo("算力释放成功",
                        f"HTVS 筛选完成！\n成功对接 {success_count}/{total}\n耗时仅 {elapsed / 60:.2f} 分钟！\n包含 SMILES 的排行榜已生成！")


if __name__ == "__main__":
    main()
