from TPP.API.top_pro_pack import Project, create_project
from pathlib import Path
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String
import pandas as pd
import json

pdb_directory = Path(r"C:\dev\rsch\clique_analysis_2021\Menv_color\Menv_color")
out_dir = Path(r"C:\dev\rsch\clique_analysis_2021\Menv_log_out\Menv_log_out")

pdb_file_names = [f.name for f in pdb_directory.iterdir() if f.suffix == ".pdb"]
out_file_names = [f.name for f in out_dir.iterdir() if f.suffix == ".out" and "id" not in f.name]
pdb_file_names.sort()
out_file_names.sort()

# print(len(pdb_file_names), len(out_file_names))


mismatch = {"pdb": [], "out": []}
data_file_names = {}
for pdb in pdb_file_names:
    pdb_clean = pdb[len("Menv_color_memb_cen_nor_"):len("Menv_color_memb_cen_nor_")+4]
    for out in out_file_names:
        out_clean = out[:4]
        if pdb_clean == out_clean:
            mismatch["pdb"].append(pdb)
            mismatch["out"].append(out)
            data_file_names[str(pdb_clean)] = (pdb, out)
            break



mismatch["pdb"].sort()
mismatch["out"].sort()

ignored_paths = [pdb for pdb in pdb_file_names if pdb not in mismatch["pdb"]]

print(len(ignored_paths), len(mismatch["pdb"]), len(pdb_file_names))
print(ignored_paths)

comment_1 = '''for out in out_file_names:
    out_clean = out[:4]
    for pdb in pdb_file_names:
        pdb_clean = pdb[len("Menv_color_memb_cen_nor_"):len("Menv_color_memb_cen_nor_")+4]
        if out_clean == pdb_clean:
            mismatch["out"].append(out)
            break
print(len(mismatch["pdb"]), len(mismatch["out"]))'''

config_path = Path.cwd() / Path("dbv5_2021_config.json")
#create_project(config_path, "dbv5_2021", pdb_directory, Path.cwd() / Path("../temp_json_dir"), ignored_paths=[pdb_directory / Path(ignore_pdb) for ignore_pdb in ignored_paths], exclude_backbone=True)
comment = '''create_project(config_path, "dbv5_2021", pdb_directory, Path.cwd() / Path("../temp_json_dir"), ignored_paths=[
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_2xok_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_3j7r_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_5t15_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_6idf_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_6idp_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_6ji8_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_6kif_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_6ilu_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_5xnl_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_6kg7_r_sc.pdb",
    "C:\\dev\\rsch\\clique_analysis_2021\\Menv_color\\Menv_color\\Menv_color_memb_cen_nor_6lqi_r_sc.pdb"

  ], exclude_backbone=True)'''
proj = Project(config_path)
proj.add_ignored_path(Path(r"C:\test_proteins\Menv_color\Menv_color\list"))
#proj.load_all_pdbs([f.stem if f not in proj.list_ignored() else "" for f in proj.list_pdb_files()])
proj.load_all_pdbs(proj.generate_default_ids())


#print(len(mismatch["pdb"]), len(proj.proteins))
#assert len(mismatch["pdb"]) == len(proj.proteins)


engine = create_engine("sqlite:///clique_traceback_v4_2021.db", echo=True)

meta = MetaData()
cliques_table = Table(
    'cliques', meta,
    Column("id", Integer, primary_key=True),
    Column("size", Integer),
    Column("clique", String),
    Column("resid", String),
    Column("oldresid", String),
    Column("layerinfo", String),
    Column("pdbname", String)
)
meta.create_all(engine)

conn = engine.connect()

def get_filtered_out_lines(out_file):
    with open(out_file, "rt") as file:
        lines = file.readlines()
        return [[i for i in line.split(" ") if i != ""] for line in lines if line.split(" ")[0].strip(" ") == "2016Menv"]


def get_layer_resid(resid, layer_ref):
    return layer_ref[resid + 1]

def get_clique_with_names_only(clique):
    clique.sort(key=lambda x: x.name)
    return ";".join([i.name for i in clique])

def get_clique_with_resid_only(clique):
    clique.sort(key=lambda x: x.name)
    return ";".join([str(i.resid) for i in clique])

def get_clique_with_old_resid_only(clique):
    clique.sort(key=lambda x: x.name)
    return ";".join([str(i.old_resid) for i in clique])

def get_clique_layer_info_only(clique, layer_ref):
    clique.sort(key=lambda x: x.name)
    resids = [res.resid for res in clique]
    return ";".join([str(get_layer_resid(resid, layer_ref)) for resid in resids])

def push_clique_to_buffer(clique, pdb_name, layer_ref, buffer):
    buffer.append({"size": len(clique), "clique": get_clique_with_names_only(clique), "resid": get_clique_with_resid_only(clique), "oldresid": get_clique_with_old_resid_only(clique),
                   "layerinfo": get_clique_layer_info_only(clique, layer_ref), "pdbname": pdb_name})

def bulk_insert_cliques_into_db(buffer, conn, table):
    return conn.execute(table.insert(), buffer)


def insert_clique_into_db(clique, pdb_name, layer_ref, conn, table):
    ins = table.insert().values(size=len(clique), clique=get_clique_with_names_only(clique),
                                resid=get_clique_with_resid_only(clique), oldresid=get_clique_with_old_resid_only(clique), layerinfo=get_clique_layer_info_only(clique, layer_ref), pdbname=pdb_name)
    result = conn.execute(ins)
    return result

bad_proteins = []


min_hydrophobic_residues = 34
residue_baseline = 30  # ALL relevant debug flags need to trigger, not just first flag thrown

for pdb_id in proj.proteins:
    pdb_id_clean = pdb_id[len("Menv_color_memb_cen_nor_"):len("Menv_color_memb_cen_nor_")+4]
    if Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean))).is_file():
        flags = [pdb_id, Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean))).__str__()]
        print("out file found for {}".format(pdb_id))
        P = proj.get_protein(pdb_id)
        hydrophobic_count = 0
        layer_ref = {}
        content = get_filtered_out_lines(Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean))))
        for line in content:
            res = line[2].strip(" ")
            id = int(line[1].strip(" "))
            layer = int(line[4].strip(" "))
            layer_ref[id] = layer
            if layer == 3 or layer == 4:
                hydrophobic_count += 1


        if hydrophobic_count < min_hydrophobic_residues:
            flags.append("below hydrophobicity baseline")
        if len(P.residues) < residue_baseline:
            flags.append("below residue baseline")
        if len(layer_ref) != len(P.residues):
            flags.append("out file / pdb residue count mismatch")

        if len(flags) > 2:
            bad_proteins.append(",".join(flags) + "\n")
        else:
            print(len(layer_ref), len(P.residues), pdb_id_clean)
            cliques = P.centroid_cliques
            buffer = []
            for clique in cliques:
                # insert_clique_into_db(clique, P.name, layer_ref, conn, cliques_table)
                push_clique_to_buffer(clique, P.name, layer_ref, buffer)
            bulk_insert_cliques_into_db(buffer, conn, cliques_table)

        # insert comment 3 here

    else:
        print("out file for {} does not exist in {}".format(pdb_id, out_dir))
        bad_proteins.append(",".join((pdb_id, "", "missing out file")) + "\n")

with open("bad_proteins_file.txt", "wt") as bp_file:
    bp_file.writelines(bad_proteins)



comment_3 = '''if hydrophobic_count >= min_hydrophobic_residues:
            P = proj.get_protein(pdb_id)
            print(len(layer_ref), len(P.residues), pdb_id_clean)
            assert len(layer_ref) == len(P.residues)
            if len(P.residues) < residue_baseline:
                print("id {} for P.name with out file {} does not meet minimum residue requirements".format(pdb_id, Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean)))))
                bad_proteins.append(",".join((pdb_id, Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean))).__str__(), "below residue baseline")) + "\n")
            else:
                cliques = P.centroid_cliques
                buffer = []
                for clique in cliques:
                    # insert_clique_into_db(clique, P.name, layer_ref, conn, cliques_table)
                    push_clique_to_buffer(clique, P.name, layer_ref, buffer)
                bulk_insert_cliques_into_db(buffer, conn, cliques_table)
        else:
            print("id {} for P.name with out file {} does not meet hydrophobicity requirements".format(pdb_id, Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean)))))
            bad_proteins.append(",".join((pdb_id, Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean))).__str__(), "below hydrophobicity baseline")) + "\n")'''