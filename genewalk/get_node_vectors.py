#!/usr/bin/env python

import os
import copy
import time
import logging
import argparse
import pickle as pkl
import networkx as nx
from multiprocessing import Pool
from genewalk.deepwalk import DeepWalk
from genewalk.nx_mg_assembler import Nx_MG_Assembler


logger = logging.getLogger('genewalk.get_node_vectors')


def run_repeat(rep, MGA):
    logger.info('%s/%s' % (rep, args.Nreps))
    DW = DeepWalk(graph=MGA)

    logger.info('generate random walks')
    start = time.time()
    DW.get_walks()
    end = time.time()
    logger.info('DW.get_walks done %.2f' % (end - start))  # in sec

    logger.info('generate node vectors')
    start = time.time()
    DW.word2vec()
    end = time.time()
    logger.info('DW.word2vec done %.2f' % (end - start))  # in sec

    # Pickle the node vectors (embeddings) and DW object
    nv = copy.deepcopy(DW.model.wv)
    filename = 'GeneWalk_DW_nv_%d.pkl' % rep
    with open(os.path.join(args.path, filename), 'wb') as f:
        pkl.dump(nv, f)
    filename = 'GeneWalk_DW_%d.pkl' % rep
    with open(os.path.join(args.path, filename), 'wb') as f:
        pkl.dump(DW, f)


if __name__ == '__main__':
    # Handle command line arguments
    parser = argparse.ArgumentParser(
        description=('Choose a path where GeneWalk files are generated '
                     '(default: ~/genewalk/ ).'))
    parser.add_argument('--path', default='~/genewalk/')
    parser.add_argument('--stmts', default='data/JQ1_HGNCidForINDRA_stmts.pkl')
    parser.add_argument('--fplx', default='data/JQ1_HGNCidForINDRA_fplx.txt')
    parser.add_argument('--path_GO', default='~/genewalk/GO/')
    parser.add_argument('--Nreps', default=10)
    args = parser.parse_args()

    # Open pickled statements
    logger.info('loading %s' % args.stmts)
    with open(os.path.join(args.path, args.stmts), 'rb') as f:
        stmts = pkl.load(f)

    logger.info('assembling network')
    MG = Nx_MG_Assembler(stmts, args.path_GO)
    del stmts

    logger.info('adding genes nodes from INDRA stmts')
    MG.MG_from_INDRA()
    MG.add_FPLXannotations(args.fplx)
    logger.info('number of gene nodes: %s' % nx.number_of_nodes(MG.graph))
    logger.info('adding GO nodes')
    MG.add_GOannotations()
    MG.add_GOontology()
    logger.info('total number of nodes in network: %s' %
                nx.number_of_nodes(MG.graph))

    # pickle the network
    filename = 'GeneWalk_MG.pkl'
    MGA = copy.deepcopy(MG.graph)
    with open(os.path.join(args.path, filename), 'wb') as f:
        pkl.dump(MGA, f, protocol=pkl.HIGHEST_PROTOCOL)
    del MG

    pool = Pool(args.nproc) if args.nproc > 1 else None
    if pool:
        run_repeat_wrapper = lambda x: run_repeat(x, MGA)
        pool.map(run_repeat_wrapper, range(1, args.Nreps + 1))
    else:
        for rep in range(1, args.Nreps + 1):
            run_repeat(rep, MG)
