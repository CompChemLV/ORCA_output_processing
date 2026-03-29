# -*- coding: utf-8 -*-
"""
Este script se encarga de recolectar datos vectoriales de cálculos químicos de ORCA.
Extrae propiedades tridimensionales / de dimensión N:
- Coordenadas XYZ de la geometría optimizada por átomo
- Cargas de Mulliken y Loewdin por átomo
- Frecuencias vibracionales
- Energías de orbitales fronteras (LUMO+3 a HOMO-3)
Genera tres archivos CSV con esquema relacional usando Molecule_Group.
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

def extraer_orbitales(ruta_out):
    """Extrae las energías de orbitales (HOMO/LUMO) desde el archivo .out."""
    orbitales = []
    if not Path(ruta_out).exists():
        return orbitales
    
    with open(ruta_out, 'r') as f:
        # readlines para procesado simple
        lines = f.readlines()
    
    in_orbital_block = False
    bloque_orbs = []
    
    for i, line in enumerate(lines):
        if "ORBITAL ENERGIES" in line:
            in_orbital_block = True
            bloque_orbs = [] # Reiniciamos si hay múltiples bloques (ej. steps de OPT), nos quedamos con el último
            continue
        
        if in_orbital_block:
            # Fin del bloque: una línea vacía, o un guion largo, o asterisco
            if ("---" in line or line.strip() == "" or "*" in line) and len(bloque_orbs) > 4:
                in_orbital_block = False
                continue
            bloque_orbs.append(line.strip())
            
    parsed_orbs = []
    for l in bloque_orbs:
        parts = l.split()
        if len(parts) >= 4 and parts[0].isdigit():
            occ = float(parts[1])
            e_ev = float(parts[3])
            parsed_orbs.append({'NO': int(parts[0]), 'OCC': occ, 'Energy_eV': e_ev})
            
    if not parsed_orbs:
        return orbitales
    
    # Encontrar el HOMO (último orbital ocupado)
    homo_idx = -1
    for i, orb in enumerate(parsed_orbs):
        if orb['OCC'] > 0:
            homo_idx = i
            
    if homo_idx != -1:
        labels = ["HOMO-3", "HOMO-2", "HOMO-1", "HOMO", "LUMO", "LUMO+1", "LUMO+2", "LUMO+3"]
        indices = [homo_idx - 3, homo_idx - 2, homo_idx - 1, homo_idx, homo_idx + 1, homo_idx + 2, homo_idx + 3, homo_idx + 4]
        
        for label, idx in zip(labels, indices):
            if 0 <= idx < len(parsed_orbs):
                orbitales.append({
                    'Orbital': label,
                    'NO': parsed_orbs[idx]['NO'],
                    'OCC': parsed_orbs[idx]['OCC'],
                    'Energy_eV': parsed_orbs[idx]['Energy_eV']
                })
    return orbitales

def cargar_propiedades_json(ruta):
    """Carga un archivo JSON de ORCA."""
    try:
        with open(ruta, 'r') as f:
            return json.load(f)
    except Exception as e:
        print("Error procesando {0}: {1}".format(ruta, e))
        return None

def agrupar_y_procesar_calculos_vectoriales(rutas):
    """
    Agrupa los cálculos y extrae:
    1. Propiedades atómicas
    2. Frecuencias
    3. Orbitales
    """
    agrupados = defaultdict(dict)
    
    # 1. Identificar archivos y clasificarlos (OPT vs FREQ vs OTRO)
    for r in rutas:
        p = Path(r)
        folder = p.parent
        base_name = p.name.replace('.property.json', '')
        
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

    datos_atomos = []
    datos_frecuencias = []
    datos_orbitales = []
    
    for key, calculos in agrupados.items():
        molecule_group = Path(key).name
        
        if "OPT" in calculos:
            main_ruta = calculos["OPT"]
        elif "SP" in calculos:
            main_ruta = calculos["SP"]
        elif "FREQ" in calculos:
            main_ruta = calculos["FREQ"]
        else:
            main_ruta = list(calculos.values())[0]
            
        data_main = cargar_propiedades_json(main_ruta)
        if not data_main: continue
        
        calc_info = data_main.get('Calculation_Info', {})
        mult = calc_info.get('Mult')
        charge = calc_info.get('Charge')

        geometrias = data_main.get('Geometries', [])
        ultima_geo = geometrias[-1] if geometrias else data_main.get('Geometry_1', {})
        if not ultima_geo and geometrias == []:
            # En caso que haya datos huérfanos sin geometries_1
            ultima_geo = data_main
            
        # --- 1. Atomos ---
        coords = ultima_geo.get('Geometry', {}).get('Coordinates', {}).get('Cartesians', [])
        mulliken_data = ultima_geo.get('Mulliken_Population_Analysis')
        loewdin_data = ultima_geo.get('Loewdin_Population_Analysis')
        
        def safely_get_charges(data):
            if not data: return []
            if isinstance(data, list) and len(data) > 0:
                return data[0].get('AtomicCharges', [])
            if isinstance(data, dict):
                return data.get('AtomicCharges', [])
            return []
            
        mulliken = safely_get_charges(mulliken_data)
        loewdin = safely_get_charges(loewdin_data)
        
        for i, coord in enumerate(coords):
            elemento = coord[0]
            if isinstance(coord[1], list):
                x, y, z = coord[1][0], coord[1][1], coord[1][2]
            else:
                x, y, z = coord[1], coord[2], coord[3]
                
            chg_m = mulliken[i] if i < len(mulliken) else None
            chg_l = loewdin[i] if i < len(loewdin) else None
            
            # Arreglo para listas de un elemento
            if isinstance(chg_m, list) and len(chg_m) > 0: chg_m = chg_m[0]
            if isinstance(chg_l, list) and len(chg_l) > 0: chg_l = chg_l[0]
            
            datos_atomos.append({
                'Molecule_Group': molecule_group,
                'Multiplicity': mult,
                'Charge': charge,
                'Atom_Index': i+1,
                'Element': elemento,
                'X': x,
                'Y': y,
                'Z': z,
                'Mulliken_Charge': chg_m,
                'Loewdin_Charge': chg_l
            })
            
        # --- 2. Frecuencias ---
        def extract_freq(json_data):
            """Extrae las frecuencias vibracionales de los datos JSON.
            Soporta formatos donde THERMOCHEMISTRY_Energies es un dict o una lista (ORCA 6.1.1).
            """
            raw_thermo = None
            if 'THERMOCHEMISTRY_Energies' in json_data:
                raw_thermo = json_data['THERMOCHEMISTRY_Energies']
            else:
                geos = json_data.get('Geometries', [])
                if geos and 'THERMOCHEMISTRY_Energies' in geos[-1]:
                    raw_thermo = geos[-1]['THERMOCHEMISTRY_Energies']
                else:
                    geom_1 = json_data.get('Geometry_1', {})
                    if 'THERMOCHEMISTRY_Energies' in geom_1:
                        raw_thermo = geom_1['THERMOCHEMISTRY_Energies']
            
            if not raw_thermo:
                return []
                
            # Si es lista, tomamos el primer elemento (ORCA 6.1.1)
            thermo_dict = raw_thermo[0] if isinstance(raw_thermo, list) and len(raw_thermo) > 0 else raw_thermo
            
            if isinstance(thermo_dict, dict):
                return thermo_dict.get('FREQ', [])
            return []
            
        freqs = extract_freq(data_main)
        if not freqs and "FREQ" in calculos and calculos["FREQ"] != main_ruta:
            d_f = cargar_propiedades_json(calculos["FREQ"])
            if d_f: freqs = extract_freq(d_f)
            
        imaginarias = []
        for i, f in enumerate(freqs):
            val = f[0] if isinstance(f, list) else f
            if val < 0:
                imaginarias.append((i+1, val))
            datos_frecuencias.append({
                'Molecule_Group': molecule_group,
                'Multiplicity': mult,
                'Charge': charge,
                'Freq_Index': i+1,
                'Frequency': val
            })
            
        if imaginarias:
            print(f">>> Aviso: Molécula '{molecule_group}' tiene frecuencias imaginarias:")
            for idx, val in imaginarias:
                print(f"    Index {idx}: {val} cm-1")
            
        # --- 3. Orbitales ---
        ruta_out = Path(main_ruta).with_suffix('').with_suffix('.out')
        if not ruta_out.exists():
            out_files = list(Path(main_ruta).parent.glob("*.out"))
            opt_outs = [o for o in out_files if 'OPT' in o.name.upper()]
            if opt_outs:
                ruta_out = opt_outs[0]
            elif out_files:
                ruta_out = out_files[0]
                
        orbs = extraer_orbitales(str(ruta_out))
        for ob in orbs:
            # Reorganizamos propiedades al diccionario
            datos_orbitales.append({
                'Molecule_Group': molecule_group,
                'Multiplicity': mult,
                'Charge': charge,
                'Orbital_Label': ob['Orbital'],
                'NO': ob['NO'],
                'OCC': ob['OCC'],
                'Energy_eV': ob['Energy_eV']
            })

    return datos_atomos, datos_frecuencias, datos_orbitales


if __name__ == "__main__":
    rutas = obtener_archivos_propiedades()
    print("Se encontraron {0} archivos de propiedades.".format(len(rutas)))

    datos_atomos, datos_frecuencias, datos_orbitales = agrupar_y_procesar_calculos_vectoriales(rutas)

    output_dir = Path("Data .csv/Propiedades Vectoriales")
    output_dir.mkdir(parents=True, exist_ok=True)

    # DataFrame de Átomos
    df_atomos = pd.DataFrame(datos_atomos)
    if not df_atomos.empty:
        nom_atomos = output_dir / "propiedades_vectoriales_atomos.csv"
        df_atomos.to_csv(nom_atomos, index=False)
        print(f"Tabla Átomos guardada: {len(df_atomos)} filas en '{nom_atomos}'")
    else:
        print("No se encontraron propiedades de Átomos válidas.")

    # DataFrame de Frecuencias
    df_freqs = pd.DataFrame(datos_frecuencias)
    if not df_freqs.empty:
        nom_freqs = output_dir / "propiedades_vectoriales_frecuencias.csv"
        df_freqs.to_csv(nom_freqs, index=False)
        print(f"Tabla Frecuencias guardada: {len(df_freqs)} filas en '{nom_freqs}'")
    else:
        print("No se encontraron Frecuencias válidas.")

    # DataFrame de Orbitales
    df_orbs = pd.DataFrame(datos_orbitales)
    if not df_orbs.empty:
        nom_orbs = output_dir / "propiedades_vectoriales_orbitales.csv"
        df_orbs.to_csv(nom_orbs, index=False)
        print(f"Tabla Orbitales guardada: {len(df_orbs)} filas en '{nom_orbs}'")
    else:
        print("No se encontraron Energías de Orbitales válidas (o no había archivos .out correspondientes).")
