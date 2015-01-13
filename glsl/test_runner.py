#!/usr/bin/env python

import os
import glob
import subprocess

passed = []
failed = []

def execTest(testfile):

    cmd = "./glslc %s" % testfile

    ret = subprocess.call(cmd, shell=True)

    if ret:
        return False
   
    return True

def stat():

    global passed
    global failed

    n = len(passed) + len(failed)

    print "# of total tests : %d" % n
    print "# of tests passed: %d" % len(passed)
    print "# of tests failed: %d" % len(failed)

    for f in failed:
        print f

def main():
    test_dir = "test"
    for f in glob.glob(os.path.join(test_dir, "*.frag")):
        ret = execTest(f)

        print "%s %s" % (f, ret)

        if ret:
            passed.append(f)
        else:
            failed.append(f)

    stat()

main()
