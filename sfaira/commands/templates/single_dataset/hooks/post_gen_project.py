import os
import shutil


def remove(filepath):
    if os.path.isfile(filepath):
        os.remove(filepath)
    elif os.path.isdir(filepath):
        shutil.rmtree(filepath)


create_extra_description = '{{ cookiecutter.create_extra_description }}' == 'True'

if not create_extra_description:
    remove('extra_description.txt')
