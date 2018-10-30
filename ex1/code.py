# # # # # # # # # # # # # # # Loading Libraries # # # # # # # # # # # # # # # # 
# Usual mathematics
import os

# Directory to build walk
path = '/home/laurits/Desktop/AU/Bachelor/3rdYear/EksFysIII/ex3'
filetype = ".txt"


for paths, dirs, files in os.walk(path):
    for name in files:
        if name.endswith(filetype):
                print(name)
