from pydirac.analysis.polarizability import get_polarizability

__all__ = ["get_atomDB"]


def get_atomDB(args):
    if isinstance(args.dir_list, str):
        dirname_lis = [args.dir_list]
    else:
        dirname_lis = args.dir_list

    for d in dirname_lis:
        dir_fullname = os.path.abspath(d)
        get_polarizability(
            dir_fullname, calc_dir_patters=args.patterns, deepth=args.deepth
        )
