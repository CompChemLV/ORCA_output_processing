# -*- coding: utf-8 -*-
"""
Este script se encarga de recolectar y consolidar datos de cálculos químicos de ORCA.
Busca archivos de propiedades (.property.json), extrae información clave como energías y tiempos,
y genera un archivo CSV resumen para facilitar el análisis posterior.
Asocia cálculos de OPT y FREQ para una misma molécula.
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

def extraer_metadata_inp(ruta_json):
    """Busca el archivo .inp correspondiente al JSON y extrae metadata de nivel de teoría."""
    p = Path(ruta_json)
    nombre_base = p.name.replace('.property.json', '')
    ruta_inp = p.with_name(f"{nombre_base}.inp")
    
    if not ruta_inp.exists():
        inps = list(p.parent.glob("*.inp"))
        if inps:
            # Si hay varios INPs y estamos en OPT, preferimos el que tenga OPT en el nombre
            inps_opt = [i for i in inps if 'OPT' in i.name.upper()]
            if inps_opt:
                ruta_inp = inps_opt[0]
            else:
                ruta_inp = inps[0]
        else:
            return {}
            
    metadata = {}
    try:
        with open(ruta_inp, 'r') as f:
            lineas = f.readlines()
            
        in_method = False
        in_basis = False
        
        for linea in lineas:
            linea = linea.strip()
            if linea.startswith('!'):
                metadata['Level_Functional'] = linea.replace('!', '').strip()
            elif linea.lower().startswith('%method'):
                in_method = True
            elif linea.lower().startswith('%basis'):
                in_basis = True
            elif linea.lower().startswith('end'):
                in_method = False
                in_basis = False
                
            if in_method:
                if linea.lower().startswith('runtyp'):
                    metadata['Level_RunTyp'] = linea.split()[1]
            if in_basis:
                if linea.lower().startswith('basis'):
                    metadata['Level_Basis'] = linea.split('"')[1] if '"' in linea else linea.split()[1]
                elif linea.lower().startswith('auxj '):
                    metadata['Level_AuxJ'] = linea.split('"')[1] if '"' in linea else linea.split()[1]
                elif linea.lower().startswith('auxjk'):
                    metadata['Level_AuxJK'] = linea.split('"')[1] if '"' in linea else linea.split()[1]
                elif linea.lower().startswith('auxc'):
                    metadata['Level_AuxC'] = linea.split('"')[1] if '"' in linea else linea.split()[1]
                    
    except Exception as e:
        print(f"Aviso: No se pudo leer {ruta_inp}: {e}")
        
    return metadata

def cargar_propiedades_json(ruta):
    """Carga un archivo JSON de ORCA y extrae el diccionario parseado."""
    try:
        with open(ruta, 'r') as f:
            return json.load(f)
    except Exception as e:
        print("Error procesando {0}: {1}".format(ruta, e))
        return None

def agrupar_y_procesar_calculos(rutas):
    """
    Agrupa los cálculos (.json) referidos a una misma molécula/carpeta,
    priorizando la data principal del cálculo de OPT (o SP) y añadiendo la 
    termoquímica del cálculo FREQ si existe.
    """
    agrupados = defaultdict(dict)
    
    # 1. Identificar archivos y clasificarlos (OPT vs FREQ vs OTRO)
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

    # 2. Procesar y consolidar cada grupo
    datos_consolidados = []
    
    for key, calculos in agrupados.items():
        # Decidir de donde tomar el perfil principal
        # Prioridad: OPT -> SP -> FREQ -> OTHER
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
        
        geometrias = data_main.get('Geometries', [])
        ultima_geo = geometrias[-1] if geometrias else {}
        geo_index = len(geometrias) - 1
        
        calc_info = data_main.get('Calculation_Info', {})
        calc_status = data_main.get('Calculation_Status', {})
        calc_timings = data_main.get('Calculation_Timings', {})
        dft_energy = ultima_geo.get('DFT_Energy', {})
        
        resumen = {
            'Molecule_Group': Path(key).name, # Sólo el nombre de la molécula
            'Main_File': Path(main_ruta).name,
            'GeometryIndex': geo_index,
            'Mult': calc_info.get('Mult'),
            'Charge': calc_info.get('Charge'),
            'NumOfAtoms': calc_info.get('NumOfAtoms'),
            'nAlphaEl': dft_energy.get('nAlphaEl'),
            'nBetaEl': dft_energy.get('nBetaEl'),
            'nTotalEl': dft_energy.get('nTotalEl'),
            'finalEn': dft_energy.get('finalEn'),
            'eExchange': dft_energy.get('eExchange'),
            'eCorr': dft_energy.get('eCorr'),
            'eXC': dft_energy.get('eXC'),
            'vdW': ultima_geo.get('VdW_Correction', {}).get('vdW'),
            'FinalEnergy': ultima_geo.get('Single_Point_Data', {}).get('FinalEnergy')
        }
        
        # Extraer Metadata de Teoría (.inp) del archivo principal de OPT/SP
        # Si el calculo que agarramos es el main_ruta, tomamos metadata de ahi
        ref_metadata_ruta = calculos.get("OPT", calculos.get("SP", main_ruta))
        metadata = extraer_metadata_inp(ref_metadata_ruta)
        for k, v in metadata.items():
            resumen[k] = v
            
        # Agregar estado y tiempos (del archivo principal)
        for k, v in calc_status.items():
            resumen[f"Status_{k}"] = v
        for k, v in calc_timings.items():
            resumen[f"Timing_{k}"] = v
            
        def extract_thermo(json_data):
            """Intenta buscar las propiedades termoquímicas tanto en la raíz como en las geometrías.
            Soporta formatos donde THERMOCHEMISTRY_Energies es un dict o una lista (ORCA 6.1.1).
            """
            raw_data = None
            if 'THERMOCHEMISTRY_Energies' in json_data:
                raw_data = json_data['THERMOCHEMISTRY_Energies']
            else:
                geos = json_data.get('Geometries', [])
                if geos and 'THERMOCHEMISTRY_Energies' in geos[-1]:
                    raw_data = geos[-1]['THERMOCHEMISTRY_Energies']
                else:
                    # Para el caso en que venga en Geometry_1
                    geom_1 = json_data.get('Geometry_1', {})
                    if 'THERMOCHEMISTRY_Energies' in geom_1:
                        raw_data = geom_1['THERMOCHEMISTRY_Energies']
            
            if isinstance(raw_data, list) and len(raw_data) > 0:
                return raw_data[0]
            return raw_data if isinstance(raw_data, dict) else None
            
        # 3. Buscar termo química: Primero intentamos extraer del archivo principal (OPT o SP que incluya FREQ)
        thermo_data = extract_thermo(data_main)
        if thermo_data:
            resumen['Freq_File'] = 'Embedded in ' + Path(main_ruta).name
        else:
            # Si no está adentro, buscamos si hay un archivo separado de FREQ
            if "FREQ" in calculos and calculos["FREQ"] != main_ruta:
                resumen['Freq_File'] = Path(calculos["FREQ"]).name
                data_freq = cargar_propiedades_json(calculos["FREQ"])
                if data_freq:
                    thermo_data = extract_thermo(data_freq)
            else:
                resumen['Freq_File'] = 'None'
            
        if thermo_data:
            resumen['Has_Freq'] = True
            for k, v in thermo_data.items():
                if not isinstance(v, (list, dict)):
                    resumen[f"Thermo_{k}"] = v
        else:
            resumen['Has_Freq'] = False

        datos_consolidados.append(resumen)
        
    return datos_consolidados


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
        print("-" * 50)
        # Mostramos las primeras 5 filas traspuestas para que se vea bien en monitores pequeños
        with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
            print(df.head(5).T)
        print("-" * 50)
        
        output_dir = Path("Data .csv/Propiedades Numericas")
        output_dir.mkdir(parents=True, exist_ok=True)
        nombre_salida = output_dir / "propiedades_numéricas_calculos.csv"
        df.to_csv(nombre_salida, index=False)
        print("\nDatos guardados exitosamente en '{0}'".format(nombre_salida))
    else:
        print("No se encontraron datos válidos para extraer.")
