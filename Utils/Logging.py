'''
Created on Jun 9, 2016
@author: lixi729

This is the class to define logging process:

print command will write to screen and a logging file
'''
import sys

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = None

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass