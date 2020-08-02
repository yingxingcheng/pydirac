import os
from pydirac.core.settings import Settings
from pydirac.input.jobs import DiracJob


def add_sub_node(settinging_obj, dir_node, subdir_node,
        keyword_node, value_list):
    """Add sub node to parent node
    """
    assert (type(dir_node) is not None)
    assert (type(value_list) == list)

    # whether add value to parent node
    # if yes: value_list has been added to parent node
    # if no: nothing has been done
    is_set = False

    if subdir_node is not None and keyword_node is not None:
        #-------------------------------------------------------
        # Case 1, for example:
        #-------------------------------------------------------
        # **ANALYZE
        # .MULPOP
        # *MULPOP
        # .VECPOP
        # 1..oo
        if len(value_list) > 0:
            settinging_obj[dir_node][subdir_node][keyword_node] = value_list
        else:
            settinging_obj[dir_node][subdir_node][keyword_node] = True
        is_set = True
    elif subdir_node is None and keyword_node is not None:
        #-------------------------------------------------------
        # Case 2, for example:
        #-------------------------------------------------------
        # **DIRAC
        # .TITLE
        # B, DOSSSS, KRCI
        if len(value_list) > 0:
            settinging_obj[dir_node][keyword_node] = value_list
        else:
            settinging_obj[dir_node][keyword_node] = True
        is_set = True
    elif subdir_node is not None and keyword_node is None:
        # if subdir_node is not None and keyword_node is None
        #-------------------------------------------------------
        # Case 3, for example:
        #-------------------------------------------------------
        # **INTEGRALS
        # *READINP
        # .UNCONTRACT
        # 
        # when it's processing .UNCONTRACT, it's the case 3
        pass
    else:
        # if subdir_node is None and keyword_node is None
        #-------------------------------------------------------
        # Case 4, for example:
        #-------------------------------------------------------
        # **DIRAC
        #
        # when processing the part in which there is only directory
        pass
    return is_set


def parse_dirac_input(file_obj='tmp.inp'):
    """
    Parse DIRAC input file to restore python Settings object
    """
    if type(file_obj) == str:
        if os.path.isfile(file_obj):
            with open(file_obj, 'r') as f:
                lines = f.readlines()
        else:
            lines = [ l + '\n' for l in file_obj.split('\n')]
    elif type(file_obj) == list:
        lines = file_obj
    else:
        raise TypeError('Type of file_obj is invalid!')

    setting = Settings()
    setting.input.dirac = Settings()

    pre_dir = None
    pre_subdir = None
    pre_dotkey = None
    curr_value_list = []

    for line in lines:
        if line.startswith('#'):
            # this is a comment line
            continue

        if line.startswith('**'):
            if pre_dir is not None:
                is_set = add_sub_node(setting.input, pre_dir, pre_subdir, pre_dotkey, curr_value_list)
                if is_set: curr_value_list = []
            if len(curr_value_list) > 0:
                add_sub_node(setting.input, pre_dir, pre_subdir, pre_dotkey, curr_value_list)

            # this is a directory
            dir_name = line.lstrip('**').rstrip()
            pre_dir = dir_name
            setting.input[pre_dir] = Settings()

            ## clear current value for subdir and keyword
            pre_subdir = None
            pre_dotkey = None
            curr_value_list.clear()

        elif line.startswith('*'):
            # this is a sub-directory
            if pre_dir is not None:
                is_set = add_sub_node(setting.input, pre_dir, pre_subdir, pre_dotkey, curr_value_list)
                if is_set: curr_value_list = []

            if line.startswith('*END'):
                continue

            sub_dir_name = line.lstrip('*').rstrip()
            pre_subdir = sub_dir_name

            assert (pre_dir is not None)
            if sub_dir_name in setting.input[pre_dir]:
                ## if sub_dir_name in pre_dir word list
                setting.input[pre_dir][sub_dir_name]  = Settings()
                setting.input[pre_dir][sub_dir_name]._en = True
            else:
                ## this is a real second level directory
                setting.input[pre_dir][pre_subdir] = Settings()

            ## clear current value for keyword
            pre_dotkey = None
            curr_value_list = []

        elif line.startswith('.'):
            # this is dot keyword
            curr_dotkey = line.lstrip('.').rstrip()
            if pre_dir is not None:
                ## first we probabily have two same dotkey, e.g., CIROOTS
                if curr_dotkey == pre_dotkey:
                    curr_value_list.append('.' + curr_dotkey)
                    continue

                is_set = add_sub_node(setting.input, pre_dir, pre_subdir, pre_dotkey, curr_value_list)
                if is_set: 
                    curr_value_list = []

                if pre_subdir is None and curr_dotkey in setting.input[pre_dir]:
                    print('{}, {}, {}'.format(pre_dir, pre_subdir, curr_dotkey))
                    curr_value_list.append('.' + curr_dotkey)
                    continue

            pre_dotkey = curr_dotkey
            assert( not (pre_dir is None and pre_subdir is None))

            ## clear value for value list
            curr_value_list = []
        else:
            # this is value of keyword
            curr_value = line.rstrip()
            if len(curr_value) == 0:
                ## empty line
                continue
            curr_value_list.append(curr_value)

    return setting


if __name__ == '__main__':
    setting = parser_dirac_input()
    job = DiracJob(settings=setting)

    with open('tmp2.inp', 'w') as f:
        f.write(job.get_input())
