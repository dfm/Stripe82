#!/usr/bin/env python
# encoding: utf-8
"""


"""

__all__ = [""]

import os
import sys
import time
import atexit
from subprocess import Popen, PIPE, STDOUT
from socket import error, socket, AF_INET, SOCK_STREAM
from threading import Thread

import pymongo
from bson.son import SON

devnull = open("/dev/null", "w+")

def waitfor(proc, port):
    tries = 0
    while proc.poll() is None and tries < 40:
        tries += 1
        s = socket(AF_INET, SOCK_STREAM)
        try:
            try:
                s.connect(('localhost', port))
                print "Started process on port %d"%port
                return
            except (IOError, error):
                time.sleep(0.25)
        finally:
            s.close()
    sys.exit(1)

def start_db(bp, port, index):
    path = os.path.join(bp, "shard_%d"%i)
    try:
        os.makedirs(path)
    except os.error:
        pass
    args = ["numactl", "--interleave=all", "mongod", "--shardsvr",
            "--dbpath", path, "--port", str(port)]
    db = Popen(args, stdin=devnull, stdout=PIPE, stderr=STDOUT)
    waitfor(db, port)
    return db

def start_config(bp, port, i, chunksize):
    path = os.path.join(bp, "config_%d"%i)
    try:
        os.makedirs(path)
    except os.error:
        pass
    args = ["numactl", "--interleave=all", "mongod", "--configsvr",
            "--dbpath", path, "--port", str(port)]
    db = Popen(args, stdin=devnull, stdout=PIPE, stderr=STDOUT)
    waitfor(db, port)
    c = pymongo.Connection("localhost", port).config
    c.settings.save({"_id": "chunksize", "value": chunksize}, safe=True)
    del c
    return db

def start_mongos(bp, port, i, configs):
    args = ["numactl", "--interleave=all", "mongos", "--configdb", configs,
            "--port", str(port)]
    db = Popen(args, stdin=devnull, stdout=PIPE, stderr=STDOUT)
    waitfor(db, port)
    return db

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--dbpath",
        help="The basepath for the databases.",
        type=str,
        default="./db/sharding")
    parser.add_argument("--ndb",
        help="The number of shards to run.",
        type=int,
        default=2)
    parser.add_argument("--nconfig",
        help="The number of config servers to run.",
        type=int,
        default=1)
    parser.add_argument("--nmongos",
        help="The number of config servers to run.",
        type=int,
        default=1)
    parser.add_argument("--chunksize",
        help="The chunk size in megabytes.",
        type=int,
        default=64)
    args = parser.parse_args()

    dbpath = os.path.abspath(args.dbpath)
    try:
        os.makedirs(dbpath)
    except os.error:
        pass

    procs = []
    fds = []
    configs = []

    @atexit.register
    def exit():
        for proc in procs:
            try:
                proc.terminate()
            except OSError:
                pass

    for i, port in enumerate(range(20000, 20000+args.nconfig)):
        config = start_config(dbpath, port, i, args.chunksize)
        config.prefix = "C%d"%i
        procs.append(config)
        configs.append("localhost:%d"%port)
    configs = ",".join(configs)

    db_ports = range(30000, 30000+args.ndb)
    for i, port in enumerate(db_ports):
        db = start_db(dbpath, port, i)
        db.prefix = "D%d"%i
        procs.append(db)

    if args.nmongos == 1:
        s_port = 27017
        ms = start_mongos(dbpath, s_port, 0, configs)
        ms.prefix = "S"
        procs.append(ms)
    else:
        s_port = 10000
        for i, port in enumerate(range(s_port, s_port+args.nmongos)):
            ms = start_mongos(dbpath, port, i, configs)
            ms.prefix = "S%d"%i
            procs.append(ms)

    # Add the shards.
    admin = pymongo.Connection("localhost", s_port).admin
    for port in db_ports:
        try:
            admin.command("addshard", "localhost:%d"%port, allowLocal=True)
        except pymongo.errors.OperationFailure:
            pass

    # Enable sharding.
    try:
        admin.command("enablesharding", "sdss")
    except pymongo.errors.OperationFailure:
        pass
    admin.command("shardcollection", "sdss.photometry",
        key=SON((k,1) for k in "run,star".split(',')))

    time.sleep(2)

    def printer():
        while len(procs):
            for proc in procs:
                line = proc.stdout.readline().rstrip()
                if line:
                    print proc.prefix, line
                else:
                    if proc.poll() is not None:
                        print proc.prefix, "EXITED", proc.returncode
                        procs.pop(proc)
                        break
                break

    printer_thread = Thread(target=printer)
    printer_thread.start()

    try:
        printer_thread.join()
    except KeyboardInterrupt:
        pass

