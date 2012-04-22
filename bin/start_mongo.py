#!/usr/bin/env python
# encoding: utf-8
"""


"""

__all__ = [""]

import os
import sys
import time
import atexit
import subprocess
from socket import error, socket, AF_INET, SOCK_STREAM

import pymongo

devnull = open("/dev/null", "w+")

def waitfor(proc, port):
    tries = 0
    while proc.poll() is None and tries < 40:
        tries += 1
        s = socket(AF_INET, SOCK_STREAM)
        try:
            try:
                s.connect(('localhost', port))
                return
            except (IOError, error):
                time.sleep(0.25)
        finally:
            s.close()
    sys.exit(1)

def start_db(bp, port, index):
    logf = open(os.path.join(bp, "shard_%d.log"%i), "a")
    path = os.path.join(bp, "shard_%d"%i)
    try:
        os.makedirs(path)
    except os.error:
        pass
    args = ["mongod", "--shardsvr", "--dbpath", path, "--port", str(port)]
    db = subprocess.Popen(args, stdin=devnull, stdout=logf, stderr=logf)
    waitfor(db, port)
    return db

def start_config(bp, port, i, chunksize):
    logf = open(os.path.join(bp, "config_%d.log"%i), "a")
    path = os.path.join(bp, "config_%d"%i)
    try:
        os.makedirs(path)
    except os.error:
        pass
    args = ["mongod", "--configsvr", "--dbpath", path, "--port", str(port)]
    db = subprocess.Popen(args, stdin=devnull, stdout=logf, stderr=logf)
    waitfor(db, port)
    c = pymongo.Connection("localhost", port).config
    c.settings.save({"_id": "chunksize", "value": chunksize}, safe=True)
    del c
    return db

def start_mongos(bp, port, i, configs):
    logf = open(os.path.join(bp, "config_%d.log"%i), "a")
    args = ["mongos", "--configdb", configs, "--port", str(port)]
    db = subprocess.Popen(args, stdin=devnull, stdout=logf, stderr=logf)
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
        procs.append(config)
        configs.append("localhost:%d"%port)
    configs = ",".join(configs)

    db_ports = range(30000, 30000+args.ndb)
    for i, port in enumerate(db_ports):
        db = start_db(dbpath, port, i)
        procs.append(db)

    if args.nmongos == 1:
        s_port = 27017
        mongos = start_mongos(dbpath, s_port, 0, configs)
        procs.append(mongos)
    else:
        s_port = 10000
        for i, port in enumerate(range(s_port, s_port+args.nmongos)):
            mongos = start_mongos(dbpath, port, i, configs)
            procs.append(mongos)

    # Add the shards.
    admin = pymongo.Connection("localhost", s_port).admin
    for port in db_ports:
        try:
            admin.command("addshard", "localhost:%d"%port, allowLocal=True)
        except pymongo.errors.OperationFailure:
            pass
    time.sleep(2)

    try:
        while True:
            pass
    except KeyboardInterrupt:
        pass

