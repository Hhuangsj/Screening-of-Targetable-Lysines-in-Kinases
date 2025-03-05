import os
import requests
import pandas as pd
import numpy as np
from pymol import cmd
from rdkit import Chem
from pymol import stored

def get_klifs_path():
    OVERVIEW_PATH = 'D:/course/yaolab/Python/klifs/overview.csv'
    FILE_PATH = 'D:/course/yaolab/Python/klifs/HUMAN'
    return OVERVIEW_PATH, FILE_PATH

def residue_around_smiles_prep(ligand_pdb_path, complex_pdb_path, query_smiles, cutoff, pymol_select):
    """
    Search a given pattern using pymol selection algebra
    Return residue name and id 
    """
    cmd.load(ligand_pdb_path, 'ligand')
    cmd.load(complex_pdb_path)
    try:
        lig = Chem.rdmolfiles.MolFromPDBFile(ligand_pdb_path)
    except:
        return 'NA'
    q = Chem.MolFromSmiles(query_smiles)

    r=[]
    try:
        match_sign = len(lig.GetSubstructMatches(q))
#        print(match_sign)
    except:
        return 'NA'
    if match_sign != 0:
        #print("pattern %s find in %s" %(query_smiles,complex_pdb_path))
        for i in range(len(lig.GetSubstructMatches(q))):
            #convert (0, 18, 17, 16, 15, 19) to 0+18+17+16+15+19
            sel_index = str(lig.GetSubstructMatches(q)[i]).replace(',','+').replace('(','').replace(')','').replace(' ','')
            cmd.select('ligAr','ligand and rank %s'%(sel_index))
            #cmd.select('lys%d'%(i),'complex and resn LYS and name NZ within %s of ligAr' %(cutoff))
            cmd.select('lys%d'%(i),'complex and %s within %s of ligAr' %(pymol_select,cutoff))
            stored.list=[]
            cmd.iterate('lys%d'%(i), 'stored.list.append((resn,resi))')
            #print(stored.list)
            try:
                r.append(stored.list[0][0]+" " +stored.list[0][1])
            except:
                a = 1
#                print("no %s within %s angstrom index '%s'" %(pymol_select,cutoff,sel_index.replace('+',' ')))
    else:
#        print("pattern %s not find in %s"%(query_smiles,complex_pdb_path) )
        r = 'NA'
    cmd.delete('all')
    return str(r).replace('[','').replace(']','').replace("'",'')

def residue_in_file_prep(file_path, **kwargs):
    """
    The ligands are prepared by preparation wizard 
    $SCHRODINGER/utilities/prepwizard ligand.pdb ligand_prep.pdb
    $SCHRODINGER/utilities/prepwizard ligand_2.pdb ligand_2_prep.pdb
    """
    if os.path.exists(file_path + "/ligand_prep.pdb"):
        ligand_pdb_path = file_path + "/ligand_prep.pdb"
    elif os.path.exists(file_path + "/ligand_2_prep.pdb"):
        ligand_pdb_path = file_path + "/ligand_2_prep.pdb"
    else:
        ligand_pdb_path = ''
#        print('ligand pdb does not exist: %s' %(file_path))

    complex_pdb_path = file_path + "/complex.pdb"
    query_smiles = kwargs['query_smiles']
    cutoff = kwargs['cutoff']
    pymol_select = kwargs['pymol_select']
    K_list = []
    if os.path.exists(ligand_pdb_path) and os.path.exists(complex_pdb_path):
        K_list = residue_around_smiles_prep(ligand_pdb_path, complex_pdb_path, query_smiles, cutoff, pymol_select)
    else:
        K_list = 'NA'
    return K_list


def lysine_around_smiles(ligand_pdb_path, complex_pdb_path, query_smiles, cutoff):
    """
    Find the lysine residue of a given PDB structure 
    The ligand PDB need to be extracted from the complex structure first
    """
    cmd.load(ligand_pdb_path)
    cmd.load(complex_pdb_path)
    lig = Chem.rdmolfiles.MolFromPDBFile(ligand_pdb_path)
    q = Chem.MolFromSmiles(query_smiles)

    r=[]
    try:
        match_sign = len(lig.GetSubstructMatches(q))
    except:
        return 'NA'
    if match_sign != 0:
        print("pattern %s find in %s" %(query_smiles,complex_pdb_path))
        for i in range(len(lig.GetSubstructMatches(q))):
            #convert (0, 18, 17, 16, 15, 19) to 0+18+17+16+15+19
            sel_index = str(lig.GetSubstructMatches(q)[i]).replace(',','+').replace('(','').replace(')','').replace(' ','')
            cmd.select('ligAr','ligand and rank %s'%(sel_index))
            cmd.select('lys%d'%(i),'complex and resn LYS and name NZ within %s of ligAr' %(cutoff))
            stored.list=[]
            cmd.iterate('lys%d'%(i), 'stored.list.append((resn,resi))')
            #print(stored.list)
            try:
                r.append(stored.list[0][0]+" " +stored.list[0][1])
            except:
                print("no lysine within %s angstrom index '%s'" %(cutoff,sel_index.replace('+',' ')))
    else:
        print("pattern %s not find in %s"%(query_smiles,complex_pdb_path) )
        r = 'NA'
    cmd.delete('all')
    return str(r).replace('[','').replace(']','').replace("'",'')



def lys_in_file(file_path, query_smiles, cutoff):
    """
    I extracted the ligand from pdb files if ligand.pdb is not exist
    """
    if os.path.exists(file_path + "/ligand.pdb"):
        ligand_pdb_path = file_path + "/ligand.pdb"
    elif os.path.exists(file_path + "/ligand_2.pdb"):
        ligand_pdb_path = file_path + "/ligand_2.pdb"
    else:
        ligand_pdb_path = ''
        print('ligand pdb does not exist: %s' %(file_path))
        
    complex_pdb_path = file_path + "/complex.pdb"
    K_list = []
    if os.path.exists(ligand_pdb_path) and os.path.exists(complex_pdb_path):
        K_list = lysine_around_smiles(ligand_pdb_path, complex_pdb_path, query_smiles, cutoff)
    else:
        K_list = 'NA'
    return K_list

def add_path_to_df(FILE_PATH,df):
    """
    Add file path to the pandas df as and return the df
    """
    file_path_list = []
    file_path = []
    base_path = FILE_PATH
    for i in range(len(df)):
        kinase = df.iloc[i]['kinase']
        pdb = df.iloc[i]['pdb']
        chain = df.iloc[i]['chain']
        alt = df.iloc[i]['alt']
        if alt == ' ':
            file_path = base_path + '/' + kinase + '/' + pdb + '_' + 'chain' + chain
        else:
            file_path = base_path + '/' + kinase + '/' + pdb + '_' + 'alt' + alt + '_' + 'chain' + chain
        file_path_list.append(file_path)

    df_file_path = df
    df_file_path['file_path'] = file_path_list
    return df_file_path

def cat_pos(unprot_id, df_catK):
    pos_lys = 'NA'
    for i in range(len(df_catK)):
        if str(unprot_id) in df_catK.iloc[i]['Accession Number  uniprot']:
            pos_lys = df_catK.iloc[i]['Position of Cat Lysine']
    return pos_lys

def check_nth_residue(pocket_pdb,n):
    """
    return the resn and resi of the nth residue in pocket.pdb 
    """
    if os.path.exists(pocket_pdb):
        cmd.load(pocket_pdb)
        myspace = {'res_list':[]}
        #iterate  CA resn and resi of pocket residues
        cmd.iterate('name ca','res_list.append(resn+" "+resi)', space=myspace)
        if len(myspace['res_list']) > 50: #without too many missing residues
            resn_resi = myspace['res_list'][n-1] #nth residue
        else:
            resn_resi = 'NA'
        cmd.delete('all')
    else:
        resn_resi = 'NA'
    return resn_resi

def correct_nth_residue_bypocket(pocket_seq,n):
    """
    count how many "_" before the nth resi
    
    6bu6 A D
    
    'LPIG__GFGCIYLCVVKVEPLFTELKFYQRAALGVPKYWGSFMIMDRFG_SDLQKIYEAYIHEHEYVHGDIKASNLLLLVDYGLA'
    
    the before the 17th residue K, there are two "_" 
    """
    count = 0
    i = 0
    while i < n:
        if pocket_seq[i] == '_':
            count += 1
        i+=1
    return count

def check_cat(lys_list, catK_byPDB, catK_byWX):
    unmatch = 0
    for i in lys_list:
        if i != catK_byPDB and i != catK_byWX:
            unmatch += 1
    if catK_byPDB == catK_byWX and unmatch != 0:
        unmatch -= 1
    return unmatch

def get_ligand_smiles(file_path):
    """
    return the SMILES of ligand.pdb
    """
    ligand_pdb_path = file_path + '/ligand.pdb'
    ligand2_pdb_path = file_path + '/ligand_2.pdb'
    if os.path.exists(ligand_pdb_path):
        lig = Chem.rdmolfiles.MolFromPDBFile(ligand_pdb_path)
        smi = Chem.MolToSmiles(lig)
    elif os.path.exists(ligand2_pdb_path):
        lig = Chem.rdmolfiles.MolFromPDBFile(ligand2_pdb_path)
        smi = Chem.MolToSmiles(lig)
    else:
        smi = 0
    return smi

def get_rcsb_info(pdbid_list):
    """
    input PDB IDs
    return json_data
    """
    url = "https://data.rcsb.org/graphql"
    
    entry_ids = ', '.join(['"' + pdbid + '"' for pdbid in pdbid_list])
    
    body = """
    {
      entries(entry_ids: [%s])
      {
    rcsb_id
    audit_author {
      name
    }
    pubmed {
      rcsb_pubmed_container_identifiers {
        pubmed_id
      }
      rcsb_pubmed_doi
    }
    rcsb_accession_info {
      initial_release_date
    }
    rcsb_entry_container_identifiers {
      entry_id
    }
    rcsb_primary_citation {
      journal_volume
      page_first
      page_last
      pdbx_database_id_DOI
      rcsb_authors
      rcsb_journal_abbrev
      title
      year
    }
    struct {
      title
    }
    struct_keywords {
      pdbx_keywords
    }
    polymer_entities {
      rcsb_entity_host_organism {
        ncbi_scientific_name
      }
      rcsb_entity_source_organism {
        ncbi_scientific_name
        ncbi_taxonomy_id
        rcsb_gene_name {
          value
        }
      }
      rcsb_polymer_entity {
        formula_weight
        rcsb_enzyme_class_combined {
          ec
        }
      }
      rcsb_polymer_entity_container_identifiers {
        reference_sequence_identifiers {
          database_accession
        }
      }
    }
    nonpolymer_entities {
      nonpolymer_comp {
        chem_comp {
          formula
          formula_weight
          id
          name
        }
        rcsb_chem_comp_descriptor {
          InChI
          InChIKey
          SMILES
        }
      }
    }
  }
}
    """ % entry_ids  # 将 % 放在字符串末尾
    # request info from rcsb
    response = requests.post(url=url, json={"query": body})
    json_data = response.json()
    return json_data

def get_dataframe_from_json(json_data):
    """
    input json_data
    return pandas dataframe
    """
    entries = json_data['data']['entries']

    # Create lists to store the data for each entry
    rcsb_ids = []
    pubmed_ids = []
    rcsb_pubmed_dois = []
    initial_release_dates = []
    journal_volumes = []
    page_firsts = []
    page_lasts = []
    rcsb_author = []
    titles = []
    years = []
    structure_titles = []
    structure_keywords = []
    ecs = []
    database_accessions = []
    gene_names = []
    journal_abbrevs = []
    audit_authors = []
    source_organisms = []
    ncbi_scientific_names = []
    taxonomy_IDs = []
    formula_weights = []
    ligand_ids = []
    ligand_names = []
    ligand_InChIs = []
    ligand_InChIKeys = []
    ligand_SMILESs = []
    ligand_formula_weights = []
    formulas = []

    # Extract relevant data from each entry
    for entry in entries:
        rcsb_id = entry['rcsb_id']
        try:
            rcsb_pubmed_doi = entry['pubmed']['rcsb_pubmed_doi']
        except:
            rcsb_pubmed_doi = 'NA'
        try:
            pubmed_id = entry['pubmed']['rcsb_pubmed_container_identifiers']['pubmed_id']
        except:
            pubmed_id = 'NA'
        #audit_author = entry['audit_author']['name']
        try:
            ncbi_scientific_name = entry['polymer_entities'][0]['rcsb_entity_host_organism'][0]['ncbi_scientific_name']
        except:
            ncbi_scientific_name = 'NA'
        try:
            source_organism = entry['polymer_entities'][0]['rcsb_entity_source_organism'][0]['ncbi_scientific_name']
        except:
            source_organism = 'NA'
        try:
            taxonomy_ID = entry['polymer_entities'][0]['rcsb_entity_source_organism'][0]['ncbi_taxonomy_id']
        except:
            taxonomy_ID = 'NA'
        structure_title = entry['struct']['title']
        structure_keyword = entry['struct_keywords']['pdbx_keywords']
        audit_author = ', '.join(item['name'] for item in entry['audit_author'])
        initial_release_date = entry['rcsb_accession_info']['initial_release_date'][0:10]
        try:
            entry['polymer_entities'][0]['rcsb_entity_source_organism'][0]['rcsb_gene_name']
            gene_name = ', '.join(item['value'] for item in entry['polymer_entities'][0]['rcsb_entity_source_organism'][0]['rcsb_gene_name'])
        except:
            gene_name = 'NA'
        try:
            journal_volume = entry['rcsb_primary_citation']['journal_volume']
        except:
            journal_volume = 'NA'
        try:
            page_first = entry['rcsb_primary_citation']['page_first']
        except:
            page_first = 'NA'
        try:
            page_last = entry['rcsb_primary_citation']['page_last']
        except:
            page_last = 'NA'
        try:
            rcsb_authors = entry['rcsb_primary_citation']['rcsb_authors']
        except:
            rcsb_authors = 'NA'
        try:
            journal_abbrev = entry['rcsb_primary_citation']['rcsb_journal_abbrev']
        except:
            journal_abbrev = 'NA'
        try:
            title = entry['rcsb_primary_citation']['title']
        except:
            title = 'NA'
        try:
            year = entry['rcsb_primary_citation']['year']
        except:
            year = 'NA'
        formula_weight = entry['polymer_entities'][0]['rcsb_polymer_entity']['formula_weight']
        try:
            entry['polymer_entities'][0]['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers'][0]['database_accession']
            database_accession = entry['polymer_entities'][0]['rcsb_polymer_entity_container_identifiers']['reference_sequence_identifiers'][0]['database_accession']
        except:
            database_accession = 'NA'
        if 'nonpolymer_entities' in entry and isinstance(entry['nonpolymer_entities'], list):
            for i in range(len(entry['nonpolymer_entities'])):
                if 'nonpolymer_comp' in entry['nonpolymer_entities'][i] and 'rcsb_chem_comp_descriptor' in entry['nonpolymer_entities'][i]['nonpolymer_comp']:
                    chem_comp_descriptor = entry['nonpolymer_entities'][i]['nonpolymer_comp']['rcsb_chem_comp_descriptor']
                    if chem_comp_descriptor is not None:
                        inchi_key = chem_comp_descriptor.get('InChIKey')
                        if isinstance(inchi_key, str):
                            ligand_InChIKey = inchi_key + ', '
                        else:
                            ligand_InChIKey = 'NA'
                    else:
                        ligand_InChIKey = 'NA'
                else:
                    ligand_InChIKey = 'NA'
        else:
            ligand_InChIKey = 'NA'
        if 'nonpolymer_entities' in entry and isinstance(entry['nonpolymer_entities'], list):
            for i in range(len(entry['nonpolymer_entities'])):
                if 'nonpolymer_comp' in entry['nonpolymer_entities'][i] and 'rcsb_chem_comp_descriptor' in entry['nonpolymer_entities'][i]['nonpolymer_comp']:
                    chem_comp_descriptor = entry['nonpolymer_entities'][i]['nonpolymer_comp']['rcsb_chem_comp_descriptor']
                    if chem_comp_descriptor is not None:
                        smiles = chem_comp_descriptor.get('SMILES')
                        if isinstance(smiles, str):
                            ligand_SMILES = smiles + '.'
                        else:
                            ligand_SMILES = 'NA'
                    else:
                        ligand_SMILES = 'NA'
                else:
                    ligand_SMILES = 'NA'
        else:
            ligand_SMILES = 'NA'
        for i in range(len(entry['nonpolymer_entities'])):
            ligand_name = entry['nonpolymer_entities'][i]['nonpolymer_comp']['chem_comp']['name'] + ', '
            ligand_id = entry['nonpolymer_entities'][i]['nonpolymer_comp']['chem_comp']['id'] + ' '
            ligand_formula_weight = str(entry['nonpolymer_entities'][i]['nonpolymer_comp']['chem_comp']['formula_weight']) + ' '
            if type(entry['nonpolymer_entities'][i]['nonpolymer_comp']['chem_comp']['formula']) is str:
                formula = entry['nonpolymer_entities'][i]['nonpolymer_comp']['chem_comp']['formula'] + ', '
            else:
                formula = 'NA'

            try:
                ligand_InChI = entry['nonpolymer_entities'][i]['nonpolymer_comp']['rcsb_chem_comp_descriptor']['InChI'] + ', '
            except:
                ligand_InChI = 'NA'
        if entry['polymer_entities'][0]['rcsb_polymer_entity']['rcsb_enzyme_class_combined']:
            ec_list = ', '.join(item['ec'] for item in entry['polymer_entities'][0]['rcsb_polymer_entity']['rcsb_enzyme_class_combined'])
        else:
            ec_list = 'NA'

        rcsb_ids.append(rcsb_id.lower())
        rcsb_pubmed_dois.append(rcsb_pubmed_doi)
        pubmed_ids.append(pubmed_id)
        initial_release_dates.append(initial_release_date)
        gene_names.append(gene_name)
        journal_volumes.append(journal_volume)
        page_firsts.append(page_first)
        page_lasts.append(page_last)
        rcsb_author.append(rcsb_authors)
        years.append(year)
        formula_weights.append(formula_weight)
        journal_abbrevs.append(journal_abbrev)
        titles.append(title)
        database_accessions.append(database_accession)
        audit_authors.append(audit_author)
        structure_titles.append(structure_title)
        structure_keywords.append(structure_keyword)
        ncbi_scientific_names.append(ncbi_scientific_name)
        source_organisms.append(source_organism)
        taxonomy_IDs.append(taxonomy_ID)
        ligand_names.append(ligand_name)
        ligand_InChIs.append(ligand_InChI)
        ligand_InChIKeys.append(ligand_InChIKey)
        ligand_SMILESs.append(ligand_SMILES)
        ligand_ids.append(ligand_id)
        formulas.append(formula)
        ligand_formula_weights.append(ligand_formula_weight)
        ecs.append(ec_list)

    # Create a DataFrame to store the data
    df_rcsb = pd.DataFrame({
        'pdb': rcsb_ids,
        'DOI':rcsb_pubmed_dois,
        'pubmed_ID':pubmed_ids,
        'Structure_Title':structure_titles,
        'Structure_Author':audit_authors,
        'Structure_Keywords':structure_keywords,
        'Accession Code(s)':database_accessions,
        'Release Date':initial_release_dates,
        'Volume':journal_volumes,
        'First_Page':page_firsts,
        'Last_Page':page_lasts,
        'Citation_Author':rcsb_author,
        'Journal_Name_Abbrev':journal_abbrevs,
        'Title':titles,
        'Year':years,
        'Molecular_Weight':formula_weights,
        'Expression_Host':ncbi_scientific_names,
        'Source_Organism':source_organisms,
        'Taxonomy_ID':taxonomy_IDs,
        'Gene_Name':gene_names,
        'Ligand_Formula':formulas,
        'Ligand_MW':ligand_formula_weights,
        'Ligand_Name':ligand_names,
        'InChI':ligand_InChIs,
        'InChIKey':ligand_InChIKeys,
        'Ligand_SMILES':ligand_SMILESs,
        'EC_Number':ecs,
        'Ligand_ID':ligand_ids
    })
    return df_rcsb

def write_format_excel(df_out, path):
    writer = pd.ExcelWriter(path, engine='xlsxwriter') 
    df_out.to_excel(writer, sheet_name='Sheet1', index=False, na_rep='NaN')

    workbook = writer.book
    worksheet = writer.sheets['Sheet1']

    for col_idx, column in enumerate(df_out.columns):
        column_length = max(df_out[column].astype(str).apply(len).max(), len(column)) + 4
        worksheet.set_column(col_idx, col_idx, column_length)

    writer.save()

def get_ligand_smiles_prep(file_path):
    """
    return the SMILES of ligand_prep.pdb
    """
    ligand_pdb_path = file_path + '/ligand_prep.pdb'
    ligand2_pdb_path = file_path + '/ligand_2_prep.pdb'
    if os.path.exists(ligand_pdb_path):
        lig = Chem.rdmolfiles.MolFromPDBFile(ligand_pdb_path)
        smi = Chem.MolToSmiles(lig)
    elif os.path.exists(ligand2_pdb_path):
        lig = Chem.rdmolfiles.MolFromPDBFile(ligand2_pdb_path)
        smi = Chem.MolToSmiles(lig)
    else:
        smi = 0
    if len(smi) > 1000:
        smi = "NA"
    return smi
    
