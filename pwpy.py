from utils.pwpy_utils import XmlPw, fort2py, py2fort

t = XmlPw()
last_tags = ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS']

def pw2py(filename):
    input_pw = {}
    for tag in last_tags:
        input_pw[tag] = {'rows': []}

    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    read = False
    for line in lines:
        sp = line.split('!')[0].split()
        if line.startswith('/'):
            read = False
        if read == 'NAMELISTS' and len(sp) == 3:
            key, _, value = sp
            if value[-1] == ',':
                value = value[:-1]
            type_of_value = t.search_by_name(key).get('type')
            input_pw[key] = fort2py(type_of_value, value)
        if line.startswith('&'):
            read = 'NAMELISTS'
        
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
                                sp = map(int, sp)
                            sp = list(map(float, sp))
                        case _:
                            raise ValueError(f'Unknown type {input_pw[read]["type"]}')
            input_pw[read]['rows'].append(sp)
    return input_pw


def py2pw(py_dict, outfile='py.pw.in'):
    for key, value in py_dict.items():
        if key not in last_tags:
            var = t.search_by_name(key)
            if not isinstance(value, XmlPw.conversions[var.get('type')]):
                raise ValueError(f'Variable {key} should be of type {var.get("type")}')
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
            lines.append(f'  {key} = {py2fort(value)}\n')
        lines.append('/\n\n')

    atomic_species = [row[0] for row in py_dict['ATOMIC_POSITIONS']['rows']]
    atomic_species = list(set(atomic_species))
    lines.append('ATOMIC_SPECIES\n')
    if 'ATOMIC_SPECIES' in py_dict:
        for row in py_dict['ATOMIC_SPECIES']['rows']:
            lines.append(' '.join(map(str, row)) + '\n')
    else:
        for species in atomic_species:
            lines.append(f'{species} {XmlPw.mass[species]} {species}.upf\n')
    lines.append('\n')

    lines.append(f'ATOMIC_POSITIONS {py_dict["ATOMIC_POSITIONS"]["type"]}\n')
    for row in py_dict['ATOMIC_POSITIONS']['rows']:
        lines.append(' '.join(map(str, row)) + '\n')
    lines.append('\n')

    lines.append(f'K_POINTS {py_dict["K_POINTS"]["type"]}\n')
    lines.append(' '.join(map(str, py_dict['K_POINTS']['rows'][0])))
    lines.append('\n')

    with open(outfile, 'w', encoding='utf-8') as f:
        f.writelines(lines)
