"""
Merge h5 data into one file.

Usage:
    merge2.py <joint_path> <set_path> [--cleanup]

Options:
    --cleanup   Delete distributed files after merging

"""

if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post

    print('0')
    args = docopt(__doc__)
    set_paths = list(pathlib.Path(args['<set_path>']).glob("*.h5"))
    #print(set_paths)
    post.merge_sets(args['<joint_path>'], set_paths, args['--cleanup'])
    print('2')
