# HPSandbox3
HPSandbox is a package of Python modules and example scripts for experimenting with the two-dimensional HP lattice model of Dill and Chan. It is ideally used as a teaching tool, or as a way to quickly prototype 2D lattice simulation ideas with easy-to-use extensible code -- a "sandbox" if you will. 

Included are HP_designing and HPSandbox. See the individual README files for more in depth descriptions. HPSandbox now has updated files that run on Python 3. The originals, written for Python 2, have been left as-is, for reference if necessary. On top of the conversion, some new features have been added to help analyze the HP sequence you are interested in. To use the new code, you should go to hp-code/HPSandbox/examples. Here, you will find the function enumerate3.py. Running "python3 enumerate3.py" or "python3 enumerate3.py -h" will bring up the help menu displaying usage and a list of the possible features. To test the code, try running "python3 enumerate3.py enumerate.conf", which will analyze the sequence PHPPHPHPPHH. For other examples, try using any other ".conf" file with enumerate3.py. To analyze a sequence you are interested in, just use the given configuration files as a template and edit as necessary.