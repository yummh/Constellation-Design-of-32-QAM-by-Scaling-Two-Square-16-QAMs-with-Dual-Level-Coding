import sys
import os
import io
import time

command_str_list = [
"simulation.exe 1.7 10.0 50 3390 215 1 >> DaS_3390.txt",
"simulation.exe 1.7 10.5 50 3390 305 1 >> DaS_3390.txt",
"simulation.exe 1.7 11.0 50 3390 295 1 >> DaS_3390.txt",
"simulation.exe 1.7 11.5 50 3390 285 1 >> DaS_3390.txt",
"simulation.exe 1.7 12.0 50 3390 245 1 >> DaS_3390.txt"
]

for command in command_str_list:
    print('\n< Start : ' + command + ' >')
    start_time = time.time()
    print(time.strftime('%Y.%m.%d - %H:%M:%S'))

    os.system(command)

    duration = time.time() - start_time
    minute = int(duration / 60)
    second = int(duration) % 60
    print("< %dminutes %dseconds >" % (minute,second))

    print('< Finish : ' + command + ' >')
    print('//////////////////////////////////////////////////////////////////////////////////////////\n')
