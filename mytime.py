# https://gist.github.com/turicas/5278558

from resource import getrusage as resource_usage, RUSAGE_SELF
from time import time as timestamp
import subprocess

if __name__ == '__main__':
    iterations = 20
    digits = 4
    
    print('run', 'real', 'sys', 'user')
    for run in range(iterations):
        start_time, start_resources = timestamp(), resource_usage(RUSAGE_SELF)
    
        p = subprocess.Popen(["./robosample", "inp.ala10"], stdout=subprocess.PIPE)
        out, errs = p.communicate()
    
        end_resources, end_time = resource_usage(RUSAGE_SELF), timestamp()

        real =  end_time - start_time
        sys = end_resources.ru_stime - start_resources.ru_stime
        user = end_resources.ru_utime - start_resources.ru_utime

        print(run + 1, round(real, digits), round(sys, digits), round(user, digits))