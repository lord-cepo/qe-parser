import os
from utils.pwpy_utils import XmlQe, fort2py, py2fort
import xml.etree.ElementTree as ET
from re import sub
last_tags = ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS', 'CELL_PARAMETERS']

def pw2py(filename):
    t = XmlQe('PW')
    input_pw = {}
    for tag in last_tags:
        input_pw[tag] = {'rows': []}

    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    read = False
    for line in lines:
        sp = line.split('!')[0]
        if line.strip().startswith('/'):
            read = False
        if read == 'NAMELISTS' and len(sp.split("=")) == 2:
            key, value = sub("[, \n]", "", sp).split("=")
            type_of_value = t.search_by_name(key).get('type')
            input_pw[key] = fort2py(type_of_value, value)
        if line.strip().startswith('&'):
            read = 'NAMELISTS'
        
        sp = sp.split()
        if any(line.startswith(tag) for tag in last_tags):
            read = sp[0]
            if len(sp) == 2:
                input_pw[read]['type'] = sp[1]
            continue
        if read in last_tags and len(sp) > 0:
            match read:
                case 'ATOMIC_SPECIES':
                    sp[1] = float(sp[1])
                case 'ATOMIC_POSITIONS':
                    sp = [sp[0]] + list(map(float, sp[1:]))
                case 'K_POINTS':
                    match input_pw[read]['type']:
                        case 'automatic':
                            sp = list(map(int, sp))
                        case 'crystal':
                            if len(sp) == 1:
                                sp = [int(sp[0])]
                            sp = list(map(float, sp))
                        case _:
                            raise ValueError(f'Unknown type {input_pw[read]["type"]}')
                case 'CELL_PARAMETERS':
                    sp = list(map(float,sp))
                case _:
                    raise NotImplementedError(f'Tag name {read} not implemented')
                
            input_pw[read]['rows'].append(sp)
    return input_pw


def py2pw(py_dict, outfile='py.pw.in', ignore_tags=[]):
    t = XmlQe('PW')
    for key, value in py_dict.items():
        if key not in last_tags:
            var = t.search_by_name(key)
            if key not in ignore_tags:
                if not isinstance(value, XmlQe.conversions[var.get('type')]):
                    raise ValueError(f"Variable {key} = {value} should be of type {XmlQe.conversions[var.get('type')]}, instead it's of type {type(value)}")
                options = [v2.strip().replace("'", '') for v in var.findall('.//opt') for v2 in v.get('val').split(',')]
                if len(options) > 0 and value not in options:
                    raise ValueError(f'Variable {value} should be one of {options}')
            
            namelist_node = var
            while namelist_node.tag != 'namelist':
                namelist_node = t.parent_map[namelist_node]
            t.namelists[namelist_node.get('name')].append((key, value))
                
    lines = []
    for namelist, variables in t.namelists.items():
        if len(variables) == 0:
            continue
        lines.append(f'&{namelist}\n')
        for key, value in variables:
            lines.append(f'  {key} = {py2fort(value) if key not in ignore_tags else value}\n')
        lines.append('/\n\n')

    lines.append('ATOMIC_SPECIES\n')
    if 'ATOMIC_SPECIES' in py_dict:
        for row in py_dict['ATOMIC_SPECIES']['rows']:
            lines.append(' '.join(map(str, row)) + '\n')
    else:
        atomic_species = list(set(atomic_species))
        atomic_species = [row[0] for row in py_dict['ATOMIC_POSITIONS']['rows']]
        for species in atomic_species:
            lines.append(f'{species} {XmlQe.mass[species]} {species}.upf\n')
    lines.append('\n')

    lines.append(f'ATOMIC_POSITIONS {py_dict["ATOMIC_POSITIONS"]["type"]}\n')
    for row in py_dict['ATOMIC_POSITIONS']['rows']:
        lines.append(' '.join(map(str, row)) + '\n')
    lines.append('\n')

    if len(py_dict['CELL_PARAMETERS']['rows']) > 0:
        lines.append(f'CELL_PARAMETERS {py_dict["CELL_PARAMETERS"]["type"]}\n')
        for row in py_dict['CELL_PARAMETERS']['rows']:
            lines.append(' '.join(map(str, row)) + '\n')
        lines.append('\n')
        
    lines.append(f'K_POINTS {py_dict["K_POINTS"]["type"]}\n')
    for kline in py_dict['K_POINTS']['rows']:
        lines.append(' '.join(map(str, kline))+'\n')
    lines.append('\n')

    with open(outfile, 'w', encoding='utf-8') as f:
        f.writelines(lines)

def out_schema(prefix, attribute):
    FILENAME = os.path.join(os.getenv("ESPRESSO_TMPDIR"), prefix+".save", "data-file-schema.xml")
    tree = ET.parse(FILENAME)
    root = tree.getroot()
    return root.find(".//"+attribute)