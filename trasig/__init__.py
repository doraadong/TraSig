import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path)

__all__ = ['main', 'trasig', 'utils']
for i in __all__:
    __import__(i)