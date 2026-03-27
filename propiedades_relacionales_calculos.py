# -*- coding: utf-8 -*-
"""
Este script se encarga de recolectar y consolidar datos de cálculos químicos de ORCA.
Busca archivos de propiedades (.property.json), extrae información clave como energías y tiempos,
y genera un archivo CSV resumen para facilitar el análisis posterior.
Asocia cálculos de OPT y FREQ para una misma molécula.
También extrae los Ordenes de Enlace de los archivos .json ("propiedades_relacionales_enlaces.csv")
y utiliza openbabel para generar las ZMatrix a partir de las coordenadas en "propiedades_vectoriales_atomos.csv".
"""
import json
import pandas as pd
from pathlib import Path
import re
from collections import defaultdict



def obtener_archivos_propiedades(directorio_base="."):
    """Busca recursivamente archivos .property.json en el directorio especificado."""
    p = Path(directorio_base)
    lista_archivos = [str(archivo.relative_to(p)) for archivo in p.rglob("*.property.json")]
    return lista_archivos



def cargar_propiedades_json(ruta):
    """Carga un archivo JSON de ORCA y extrae el diccionario parseado."""
    try:
        with open(ruta, 'r') as f:
            return json.load(f)
    except Exception as e:
        print("Error procesando {0}: {1}".format(ruta, e))
        return None


def agrupar_rutas_por_tipo(rutas):
    """
    Agrupa los archivos .property.json por molécula/carpeta y clasifica por tipo (OPT, FREQ, SP).
    """
    agrupados = defaultdict(dict)
    for r in rutas:
        p = Path(r)
        folder = p.parent
        base_name = p.name.replace('.property.json', '')
        
        # Identificar tipo por sufijo en nombre
        tipo = "OTHER"
        if "_OPT" in base_name.upper():
            mol_name = re.sub(r'(_OPT|_opt).*', '', base_name)
            tipo = "OPT"
        elif "_FREQ" in base_name.upper():
            mol_name = re.sub(r'(_FREQ|_freq).*', '', base_name)
            tipo = "FREQ"
        elif "_SP" in base_name.upper():
            mol_name = re.sub(r'(_SP|_sp).*', '', base_name)
            tipo = "SP"
        else:
            mol_name = base_name
        
        key = str(folder / mol_name)
        agrupados[key][tipo] = r
    return agrupados

def obtener_archivo_principal(calculos):
    """
    Decide de donde tomar el perfil principal.
    Prioridad: OPT -> SP -> FREQ -> OTHER
    """
    if "OPT" in calculos:
        return calculos["OPT"], "OPT"
    elif "SP" in calculos:
        return calculos["SP"], "SP"
    elif "FREQ" in calculos:
        return calculos["FREQ"], "FREQ"
    elif calculos:
        tipo = list(calculos.keys())[0]
        return calculos[tipo], tipo
    return None, None

def agrupar_y_procesar_calculos(rutas):
    """
    Agrupa los cálculos (.json) referidos a una misma molécula/carpeta,
    priorizando la data principal del cálculo de OPT (o SP) y añadiendo la 
    termoquímica del cálculo FREQ si existe.
    """
    agrupados = agrupar_rutas_por_tipo(rutas)

    # 2. Procesar y consolidar cada grupo
    datos_consolidados = []
    
    for key, calculos in agrupados.items():
        main_ruta, _ = obtener_archivo_principal(calculos)
        if not main_ruta: continue
            
        data_main = cargar_propiedades_json(main_ruta)
        if not data_main: continue
        
        calc_info = data_main.get('Calculation_Info', {})
        mult = calc_info.get('Mult')
        charge = calc_info.get('Charge')

        geometrias = data_main.get('Geometries', [])
        ultima_geo = geometrias[-1] if geometrias else data_main.get('Geometry_1', {})
        geo_index = len(geometrias) - 1 if geometrias else 0
        
        molecule_group = Path(key).name
        
        resumen = {
            'Molecule_Group': molecule_group,
            'Multiplicity': mult,
            'Charge': charge,
            'Main_File': Path(main_ruta).name,
            'GeometryIndex': geo_index
        }

        datos_consolidados.append(resumen)
        

    return datos_consolidados

def extraer_coordenadas_internas_outs(rutas_json):
    """
    Lee las coordenadas internas (Z-Matrix) guardadas en la salida .out de ORCA.
    """
    agrupados = agrupar_rutas_por_tipo(rutas_json)

    todos_zmat = []
    
    for key, calculos in agrupados.items():
        main_ruta, tipo_main = obtener_archivo_principal(calculos)
        if not main_ruta: continue
            
        folder = Path(main_ruta).parent
        ruta_out = None
        
        for f in folder.glob("*.out"):
            if f"_{tipo_main}.out" in f.name or f"_{tipo_main}.OUT" in f.name.upper():
                ruta_out = f
                break
                
        if not ruta_out:
            outs = list(folder.glob("*.out"))
            if outs:
                ruta_out = outs[0]
        
        if not ruta_out or not ruta_out.exists():
            print(f"No se encontró archivo .out en {folder}")
            continue
            
        molecule_group = Path(key).name
        
        try:
            with open(ruta_out, 'r', encoding='utf-8') as f:
                lines = f.readlines()
        except Exception as e:
            print(f"Error leyendo {ruta_out}: {e}")
            continue
            
        # Encontrar la última coincidencia de "INTERNAL COORDINATES (ANGSTROEM)"
        start_idx = -1
        for i in range(len(lines)-1, -1, -1):
            if "INTERNAL COORDINATES (ANGSTROEM)" in lines[i]:
                start_idx = i + 2
                break
                
        if start_idx == -1:
            continue
            
        mapping_elementos = {}
        # Primero recolectamos los elementos para poder asociarlos a las distancias/ángulos
        for i in range(start_idx, len(lines)):
            line = lines[i].strip()
            if not line or "---" in line:
                break
            parts = line.split()
            if len(parts) >= 7:
                atom_idx = i - start_idx + 1
                mapping_elementos[atom_idx] = parts[0]

        # Cargar JSON para obtener mult y charge
        data_json = cargar_propiedades_json(main_ruta)
        mult = data_json.get('Calculation_Info', {}).get('Mult') if data_json else None
        charge = data_json.get('Calculation_Info', {}).get('Charge') if data_json else None

        # Luego extraemos toda la info
        for i in range(start_idx, len(lines)):
            line = lines[i].strip()
            if not line or "---" in line:
                break
            parts = line.split()
            if len(parts) >= 7:
                element = parts[0]
                idx_bond = int(parts[1])
                idx_angl = int(parts[2])
                idx_dihe = int(parts[3])
                bond_length = float(parts[4])
                bond_angle = float(parts[5])
                dihe_angle = float(parts[6])
                
                atom_idx = i - start_idx + 1
                
                todos_zmat.append({
                    'Molecule_Group': molecule_group,
                    'Multiplicity': mult,
                    'Charge': charge,
                    'Atom_Index': atom_idx,
                    'Element': element,
                    'Bond_To': idx_bond,
                    'Bond_To_Element': mapping_elementos.get(idx_bond, ""),
                    'Angle_To': idx_angl,
                    'Angle_To_Element': mapping_elementos.get(idx_angl, ""),
                    'Dihedral_To': idx_dihe,
                    'Dihedral_To_Element': mapping_elementos.get(idx_dihe, ""),
                    'Bond_Length': bond_length,
                    'Bond_Angle': bond_angle,
                    'Dihedral_Angle': dihe_angle
                })

    return todos_zmat

def guardar_csvs_internos(datos_zmatrix):
    if not datos_zmatrix:
        print("No se extrajeron coordenadas internas válidas.")
        return
        
    df_zmat = pd.DataFrame(datos_zmatrix)
    
    output_dir = Path("Data .csv/Propiedades Relacionales")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Guardar completo
    nombre_zmat = output_dir / "propiedades_relacionales_coordenadas_internas.csv"
    df_zmat.to_csv(nombre_zmat, index=False)
    print(f"Tabla Coordenadas Internas guardada: {len(df_zmat)} filas en '{nombre_zmat}'")
    
    # 1. Distancias
    df_dist = df_zmat[df_zmat['Bond_To'] > 0][['Molecule_Group', 'Atom_Index', 'Element', 'Bond_To', 'Bond_To_Element', 'Bond_Length']].copy()
    df_dist.to_csv(output_dir / "propiedades_relacionales_distancias_enlace.csv", index=False)
    
    # 2. Angulos
    df_ang = df_zmat[df_zmat['Angle_To'] > 0][['Molecule_Group', 'Atom_Index', 'Element', 'Bond_To', 'Bond_To_Element', 'Angle_To', 'Angle_To_Element', 'Bond_Angle']].copy()
    df_ang.to_csv(output_dir / "propiedades_relacionales_angulos_enlace.csv", index=False)
    
    # 3. Diedros
    df_dihe = df_zmat[df_zmat['Dihedral_To'] > 0][['Molecule_Group', 'Atom_Index', 'Element', 'Bond_To', 'Bond_To_Element', 'Angle_To', 'Angle_To_Element', 'Dihedral_To', 'Dihedral_To_Element', 'Dihedral_Angle']].copy()
    df_dihe.to_csv(output_dir / "propiedades_relacionales_angulos_diedros.csv", index=False)
    
    print("Tablas separadas para distancias, ángulos y diedros guardadas con éxito.")

if __name__ == "__main__":
    rutas = obtener_archivos_propiedades()
    print("Se encontraron {0} archivos de propiedades.".format(len(rutas)))

    datos_extraidos = agrupar_y_procesar_calculos(rutas)

    df = pd.DataFrame(datos_extraidos)
    
    if not df.empty:
        cols = df.columns.tolist()
        status_cols = sorted([c for c in cols if c.startswith('Status_')])
        timing_cols = sorted([c for c in cols if c.startswith('Timing_')])
        thermo_cols = sorted([c for c in cols if c.startswith('Thermo_')])
        level_cols = sorted([c for c in cols if c.startswith('Level_')])
        resto_cols = [c for c in cols if not c.startswith('Status_') and not c.startswith('Timing_') and not c.startswith('Thermo_') and not c.startswith('Level_')]
        
        nuevo_orden = status_cols + level_cols + timing_cols + thermo_cols + resto_cols
        df = df[nuevo_orden]

        print(f"\nResumen de datos agrupados ({len(df)} moléculas distintas):")
        with pd.option_context('display.max_columns', None, 'display.width', 1000):
            print(df.head().to_string(index=False))
        
        output_dir_main = Path("Data .csv/Propiedades Relacionales")
        output_dir_main.mkdir(parents=True, exist_ok=True)
        nombre_salida = output_dir_main / "propiedades_relacionales_calculos.csv"
        df.to_csv(nombre_salida, index=False)
        print(f"\nDatos relacionales guardados exitosamente en '{nombre_salida}'")
    else:
        print("No se encontraron datos válidos para extraer propiedades relacionales generales.")
        

        
    # Extraer y guardar Coordenadas Internas de las salidas .out
    datos_zmatrix = extraer_coordenadas_internas_outs(rutas)
    guardar_csvs_internos(datos_zmatrix)
