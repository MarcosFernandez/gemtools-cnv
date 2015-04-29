#!/usr/bin/env python
"""Gemtools commands"""
import argparse
import gem
import gem.utils
import gem.production
import sys
from sys import exit
import signal

__VERSION__ = "1.7.1.0.0"


def gemtools():
    try:
        parser = argparse.ArgumentParser(prog="gemtools-cnv",
                description="Gemtools-cnv driver to execute different gemtools-cnv command and pipelines"
                )
        parser.add_argument('--loglevel', dest="loglevel", default=None, help="Log level (error, warn, info, debug)")
        parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __VERSION__)

        commands = {
            "index": gem.production.Index,
            "rna-pipeline": gem.production.RnaPipeline,
            "basic-mapping": gem.production.BasicMappingPipeline,
            "remove-duplicates": gem.production.RemoveDuplicatesPipeline,            
            "bam-2-fastq": gem.production.BamToFastqPipeline,
            "fastq-splitter": gem.production.SplitFastqPipeline,
            "cnv-mapping": gem.production.CnvMappingPipeline,
            "cnv-call": gem.production.CnvCallingPipeline,
            "cnv-assembly-preparation": gem.production.CnvAssemblyPreparation,
            "cnv-prep": gem.production.CnvPrepare,
            "t-index": gem.production.TranscriptIndex,
            "gtfcount": gem.production.GtfCount,
            "merge": gem.production.Merge,
            "convert": gem.production.Convert,
            "gtf-junctions": gem.production.Junctions,
            "denovo-junctions": gem.production.JunctionExtraction,
            "stats": gem.production.Stats,
            "filter": gem.production.Filter,
            "report": gem.production.StatsReport,
        }
        instances = {}

        
        subparsers = parser.add_subparsers(title="commands", metavar="<command>", description="Available commands", dest="command")
        
        for name, cmdClass in commands.items():
            p = subparsers.add_parser(name, help=cmdClass.title, description=cmdClass.description)
            instances[name] = cmdClass()
            instances[name].register(p)

        args = parser.parse_args()
        if args.loglevel is not None:
            gem.loglevel(args.loglevel)
            

        try:
            instances[args.command].run(args)
        except gem.utils.CommandException, e:
            sys.stderr.write("%s\n" % (str(e)))
            exit(1)
    except KeyboardInterrupt:
        exit(1)
    finally:
        pass

if __name__ == "__main__":
    gemtools()
