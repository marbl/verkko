MODULE       := verkko
VERSION      := snapshot v2.0
VERSION_H    := ${CANUROOT}utility/src/version.H

TARGET       := libverkko.a
SOURCES      := ${CANUROOT}utility/src/align/align-ksw2-driver.C \
                ${CANUROOT}utility/src/align/align-ksw2-extz.C \
                ${CANUROOT}utility/src/align/align-ksw2-extz2-sse.C \
                ${CANUROOT}utility/src/align/align-parasail-driver.C \
                ${CANUROOT}utility/src/align/align-ssw-driver.C \
                ${CANUROOT}utility/src/align/align-ssw.C \
                ${CANUROOT}utility/src/align/edlib.C \
                \
                ${CANUROOT}utility/src/bits/fibonacci-v1.C \
                ${CANUROOT}utility/src/bits/hexDump-v1.C \
                ${CANUROOT}utility/src/bits/stuffedBits-v1-binary.C \
                ${CANUROOT}utility/src/bits/stuffedBits-v1-bits.C \
                ${CANUROOT}utility/src/bits/stuffedBits-v1-delta.C \
                ${CANUROOT}utility/src/bits/stuffedBits-v1-gamma.C \
                ${CANUROOT}utility/src/bits/stuffedBits-v1-golomb.C \
                ${CANUROOT}utility/src/bits/stuffedBits-v1-omega.C \
                ${CANUROOT}utility/src/bits/stuffedBits-v1-unary.C \
                ${CANUROOT}utility/src/bits/stuffedBits-v1-zeckendorf.C \
                ${CANUROOT}utility/src/bits/stuffedBits-v1.C \
                ${CANUROOT}utility/src/bits/wordArray-v1.C \
                \
                ${CANUROOT}utility/src/datastructures/keyAndValue-v1.C \
                ${CANUROOT}utility/src/datastructures/splitToWords-v1.C \
                ${CANUROOT}utility/src/datastructures/stringList-v1.C \
                ${CANUROOT}utility/src/datastructures/strings-v1.C \
                ${CANUROOT}utility/src/datastructures/types-v1.C \
                \
                ${CANUROOT}utility/src/files/accessing-v1.C \
                ${CANUROOT}utility/src/files/buffered-v1-reading.C \
                ${CANUROOT}utility/src/files/buffered-v1-writing.C \
                ${CANUROOT}utility/src/files/compressed-v1-reading.C \
                ${CANUROOT}utility/src/files/compressed-v1-writing.C \
                ${CANUROOT}utility/src/files/compressed-v1.C \
                ${CANUROOT}utility/src/files/fasta-fastq-v1.C \
                ${CANUROOT}utility/src/files/files-v1.C \
                ${CANUROOT}utility/src/files/memoryMapped-v1.C \
                ${CANUROOT}utility/src/files/readLine-v0.C \
                ${CANUROOT}utility/src/files/readLine-v1.C \
                ${CANUROOT}utility/src/files/reading-v1.C \
                ${CANUROOT}utility/src/files/writing-v1.C \
                \
                ${CANUROOT}utility/src/kmers-v1/kmers-exact.C \
                ${CANUROOT}utility/src/kmers-v1/kmers-files.C \
                ${CANUROOT}utility/src/kmers-v1/kmers-histogram.C \
                ${CANUROOT}utility/src/kmers-v1/kmers-reader.C \
                ${CANUROOT}utility/src/kmers-v1/kmers-writer-block.C \
                ${CANUROOT}utility/src/kmers-v1/kmers-writer-stream.C \
                ${CANUROOT}utility/src/kmers-v1/kmers-writer.C \
                ${CANUROOT}utility/src/kmers-v1/kmers.C \
                \
                ${CANUROOT}utility/src/kmers-v2/kmers-exact.C \
                ${CANUROOT}utility/src/kmers-v2/kmers-files.C \
                ${CANUROOT}utility/src/kmers-v2/kmers-histogram.C \
                ${CANUROOT}utility/src/kmers-v2/kmers-reader-dump.C \
                ${CANUROOT}utility/src/kmers-v2/kmers-reader.C \
                ${CANUROOT}utility/src/kmers-v2/kmers-writer-block.C \
                ${CANUROOT}utility/src/kmers-v2/kmers-writer-stream.C \
                ${CANUROOT}utility/src/kmers-v2/kmers-writer.C \
                ${CANUROOT}utility/src/kmers-v2/kmers.C \
                \
                ${CANUROOT}utility/src/math/md5-v1.C \
                ${CANUROOT}utility/src/math/mt19937ar-v1.C \
                ${CANUROOT}utility/src/math/sampledDistribution-v1.C \
                \
                ${CANUROOT}utility/src/parasail/cpuid.c \
                ${CANUROOT}utility/src/parasail/memory.c \
                ${CANUROOT}utility/src/parasail/sg.c \
                ${CANUROOT}utility/src/parasail/sg_trace.c \
                ${CANUROOT}utility/src/parasail/sg_qx_dispatch.c \
                ${CANUROOT}utility/src/parasail/sg_qb_de_dispatch.c \
                ${CANUROOT}utility/src/parasail/sg_qe_db_dispatch.c \
                ${CANUROOT}utility/src/parasail/cigar.c \
                \
                ${CANUROOT}utility/src/sequence/dnaSeq-v1.C \
                ${CANUROOT}utility/src/sequence/dnaSeqFile-v1.C \
                ${CANUROOT}utility/src/sequence/sequence-v1.C \
                \
                ${CANUROOT}utility/src/system/logging-v1.C \
                ${CANUROOT}utility/src/system/runtime-v1.C \
                ${CANUROOT}utility/src/system/speedCounter-v1.C \
                ${CANUROOT}utility/src/system/sweatShop-v1.C \
                ${CANUROOT}utility/src/system/system-stackTrace-v1.C \
                ${CANUROOT}utility/src/system/system-v1.C \
                ${CANUROOT}utility/src/system/time-v1.C \
                \
                ${CANUROOT}stores/sqCache.C \
                ${CANUROOT}stores/sqLibrary.C \
                ${CANUROOT}stores/sqReadData.C \
                ${CANUROOT}stores/sqReadDataWriter.C \
                ${CANUROOT}stores/sqStore.C \
                ${CANUROOT}stores/sqStoreBlob.C \
                ${CANUROOT}stores/sqStoreConstructor.C \
                ${CANUROOT}stores/sqStoreInfo.C \
                \
                ${CANUROOT}stores/ovOverlap.C \
                ${CANUROOT}stores/ovStore.C \
                ${CANUROOT}stores/ovStoreWriter.C \
                ${CANUROOT}stores/ovStoreFilter.C \
                ${CANUROOT}stores/ovStoreFile.C \
                ${CANUROOT}stores/ovStoreHistogram.C \
                \
                ${CANUROOT}stores/tgStore.C \
                ${CANUROOT}stores/tgTig.C \
                ${CANUROOT}stores/tgTigSizeAnalysis.C \
                ${CANUROOT}stores/tgTigMultiAlignDisplay.C \
                \
                ${CANUROOT}stores/libsnappy/snappy-sinksource.cc \
                ${CANUROOT}stores/libsnappy/snappy-stubs-internal.cc \
                ${CANUROOT}stores/libsnappy/snappy.cc \
                \
                ${CANUROOT}stores/objectStore.C \
                \
                ${CANUROOT}overlapInCore/liboverlap/Binomial_Bound.C \
                \
                ${CANUROOT}gfa/gfa.C \
                ${CANUROOT}gfa/bed.C


ifeq (${BUILDSTACKTRACE}, 1)
SOURCES      += ${CANUROOT}utility/src/system/libbacktrace/atomic.c \
                ${CANUROOT}utility/src/system/libbacktrace/backtrace.c \
                ${CANUROOT}utility/src/system/libbacktrace/dwarf.c \
                ${CANUROOT}utility/src/system/libbacktrace/elf.c \
                ${CANUROOT}utility/src/system/libbacktrace/fileline.c \
                ${CANUROOT}utility/src/system/libbacktrace/mmap.c \
                ${CANUROOT}utility/src/system/libbacktrace/mmapio.c \
                ${CANUROOT}utility/src/system/libbacktrace/posix.c \
                ${CANUROOT}utility/src/system/libbacktrace/print.c \
                ${CANUROOT}utility/src/system/libbacktrace/simple.c \
                ${CANUROOT}utility/src/system/libbacktrace/sort.c \
                ${CANUROOT}utility/src/system/libbacktrace/state.c \
                ${CANUROOT}utility/src/system/libbacktrace/unknown.c
endif



SRC_INCDIRS  := . \
                ${CANUROOT}utility/src \
                ${CANUROOT}stores \
                ${CANUROOT}stores/libsnappy \
                ${CANUROOT}alignment \
                ${CANUROOT}utgcns/libcns \
                ${CANUROOT}utgcns/libpbutgcns \
                ${CANUROOT}overlapBasedTrimming \

SUBMAKEFILES := ${CANUROOT}stores/ovStoreBuild.mk \
                ${CANUROOT}stores/ovStoreConfig.mk \
                ${CANUROOT}stores/ovStoreBucketizer.mk \
                ${CANUROOT}stores/ovStoreSorter.mk \
                ${CANUROOT}stores/ovStoreIndexer.mk \
                ${CANUROOT}stores/sqStoreCreate.mk \
                ${CANUROOT}stores/sqStoreDumpMetaData.mk \
                ${CANUROOT}stores/sqStoreDumpFASTQ.mk \
                \
                ${CANUROOT}meryl/src/meryl/meryl.mk \
                ${CANUROOT}meryl/src/meryl-lookup/meryl-lookup.mk \
                \
                ${CANUROOT}overlapInCore/overlapImport.mk \
                \
                ${CANUROOT}overlapErrorAdjustment/findErrors.mk \
                ${CANUROOT}overlapErrorAdjustment/fixErrors.mk \
                \
                ${CANUROOT}utgcns/layoutToPackage.mk \
                ${CANUROOT}utgcns/utgcns.mk \
                \
                ${CANUROOT}gfa/alignGFA.mk

EXECUTABLES  := verkko.sh                             -> ../bin/verkko \
                \
                ${TARGET_DIR}/bin/findErrors          -> ../lib/verkko/bin/findErrors \
                ${TARGET_DIR}/bin/fixErrors           -> ../lib/verkko/bin/fixErrors \
                ${TARGET_DIR}/bin/layoutToPackage     -> ../lib/verkko/bin/layoutToPackage \
                ${TARGET_DIR}/bin/meryl               -> ../lib/verkko/bin/meryl \
                ${TARGET_DIR}/bin/meryl-lookup        -> ../lib/verkko/bin/meryl-lookup \
                ${TARGET_DIR}/bin/ovStoreBuild        -> ../lib/verkko/bin/ovStoreBuild \
                ${TARGET_DIR}/bin/ovStoreConfig       -> ../lib/verkko/bin/ovStoreConfig \
                ${TARGET_DIR}/bin/overlapImport       -> ../lib/verkko/bin/overlapImport \
                ${TARGET_DIR}/bin/sqStoreCreate       -> ../lib/verkko/bin/sqStoreCreate \
                ${TARGET_DIR}/bin/sqStoreDumpMetaData -> ../lib/verkko/bin/sqStoreDumpMetaData \
                ${TARGET_DIR}/bin/sqStoreDumpFASTQ    -> ../lib/verkko/bin/sqStoreDumpFASTQ \
                ${TARGET_DIR}/bin/utgcns              -> ../lib/verkko/bin/utgcns \
                ${TARGET_DIR}/bin/alignGFA            -> ../lib/verkko/bin/alignGFA

FILES        += \
                data/human-ebv-AJ507799.2.fasta.gz            -> ../lib/verkko/data/human-ebv-AJ507799.2.fasta.gz \
                data/human-mito-NC_012920.1.fasta.gz          -> ../lib/verkko/data/human-mito-NC_012920.1.fasta.gz \
                data/human-rdna-KY962518.1.fasta.gz           -> ../lib/verkko/data/human-rdna-KY962518.1.fasta.gz \
                data/chm13_rDNAs.fa                           -> ../lib/verkko/data/chm13_rDNAs.fa \
                \
                scripts/add_fake_alignments.py                -> ../lib/verkko/scripts/add_fake_alignments.py \
                scripts/add_fake_bridging_paths.py            -> ../lib/verkko/scripts/add_fake_bridging_paths.py \
                scripts/add_hairpin_uniques.py                -> ../lib/verkko/scripts/add_hairpin_uniques.py \
                scripts/add_hom_node_scaffold_edges.py        -> ../lib/verkko/scripts/add_hom_node_scaffold_edges.py \
                scripts/find_tips.py                          -> ../lib/verkko/scripts/find_tips.py \
                scripts/unroll_simple_loops.py                -> ../lib/verkko/scripts/unroll_simple_loops.py \
                scripts/select_best_alignment.py              -> ../lib/verkko/scripts/select_best_alignment.py \
                scripts/select_unused_reads.py                -> ../lib/verkko/scripts/select_unused_reads.py \
                scripts/filter_alignments_by_column.py        -> ../lib/verkko/scripts/filter_alignments_by_column.py \
                scripts/calculate_coverage.py                 -> ../lib/verkko/scripts/calculate_coverage.py \
                scripts/check_layout_gaps.py                  -> ../lib/verkko/scripts/check_layout_gaps.py \
                scripts/chop_misassemblies.py                 -> ../lib/verkko/scripts/chop_misassemblies.py \
                scripts/connect_uniques.py                    -> ../lib/verkko/scripts/connect_uniques.py \
                scripts/estimate_unique_local.py              -> ../lib/verkko/scripts/estimate_unique_local.py \
                scripts/fasta_combine.py                      -> ../lib/verkko/scripts/fasta_combine.py \
                scripts/fasta_extract.py                      -> ../lib/verkko/scripts/fasta_extract.py \
                scripts/fasta_filter.py                       -> ../lib/verkko/scripts/fasta_filter.py \
                scripts/fasta_partition.py                    -> ../lib/verkko/scripts/fasta_partition.py \
                scripts/fasta_util.py                         -> ../lib/verkko/scripts/fasta_util.py \
                scripts/find_bridges.py                       -> ../lib/verkko/scripts/find_bridges.py \
                scripts/fix_diploid_paths.py                  -> ../lib/verkko/scripts/fix_diploid_paths.py \
                scripts/fix_diploid_unique_nodes.py           -> ../lib/verkko/scripts/fix_diploid_unique_nodes.py \
                scripts/fix_haplogaps.py                      -> ../lib/verkko/scripts/fix_haplogaps.py \
                scripts/forbid_unbridged_tangles.py           -> ../lib/verkko/scripts/forbid_unbridged_tangles.py \
                scripts/get_bridge_mapping.py                 -> ../lib/verkko/scripts/get_bridge_mapping.py \
                scripts/get_layout_from_mbg.py                -> ../lib/verkko/scripts/get_layout_from_mbg.py \
                scripts/get_original_coverage.py              -> ../lib/verkko/scripts/get_original_coverage.py \
                scripts/get_unroll_mapping.py                 -> ../lib/verkko/scripts/get_unroll_mapping.py \
                scripts/inject_coverage.py                    -> ../lib/verkko/scripts/inject_coverage.py \
                scripts/insert_aln_gaps.py                    -> ../lib/verkko/scripts/insert_aln_gaps.py \
                scripts/maybe_trim_alignment.py               -> ../lib/verkko/scripts/maybe_trim_alignment.py \
                scripts/merge_unresolved_dbg_nodes.py         -> ../lib/verkko/scripts/merge_unresolved_dbg_nodes.py \
                scripts/pick_majority_bridge.py               -> ../lib/verkko/scripts/pick_majority_bridge.py \
                scripts/pop_bubbles_coverage_based.py         -> ../lib/verkko/scripts/pop_bubbles_coverage_based.py \
                scripts/pop_bubbles_keep_longest.py           -> ../lib/verkko/scripts/pop_bubbles_keep_longest.py \
                scripts/remove_contained_spurious_uniques.py  -> ../lib/verkko/scripts/remove_contained_spurious_uniques.py \
                scripts/remove_crosslink_paths.py             -> ../lib/verkko/scripts/remove_crosslink_paths.py \
                scripts/remove_wrong_connections_2.py         -> ../lib/verkko/scripts/remove_wrong_connections_2.py \
                scripts/resolve_triplets_kmerify.py           -> ../lib/verkko/scripts/resolve_triplets_kmerify.py \
                scripts/screen-assembly.pl                    -> ../lib/verkko/scripts/screen-assembly.pl \
                scripts/trim_dbg_alignment.py                 -> ../lib/verkko/scripts/trim_dbg_alignment.py \
                scripts/unitigify.py                          -> ../lib/verkko/scripts/unitigify.py \
                scripts/unroll_tip_loops.py                   -> ../lib/verkko/scripts/unroll_tip_loops.py \
                scripts/untip_relative.py                     -> ../lib/verkko/scripts/untip_relative.py \
                scripts/circularize_ctgs.py                   -> ../lib/verkko/scripts/circularize_ctgs.py \
                scripts/cluster.py                            -> ../lib/verkko/scripts/cluster.py \
                scripts/graph_functions.py                    -> ../lib/verkko/scripts/graph_functions.py \
                scripts/hicverkko.py                          -> ../lib/verkko/scripts/hicverkko.py \
                scripts/parse_sam_pairs.py                    -> ../lib/verkko/scripts/parse_sam_pairs.py \
                scripts/rdna_scaff_functions.py               -> ../lib/verkko/scripts/rdna_scaff_functions.py \
                scripts/rdna_scaff.py                         -> ../lib/verkko/scripts/rdna_scaff.py \
                scripts/remove_nodes_add_telomere.py          -> ../lib/verkko/scripts/remove_nodes_add_telomere.py \
                scripts/scaffolding/logger_wrap.py            -> ../lib/verkko/scripts/scaffolding/logger_wrap.py \
                scripts/scaffolding/match_graph.py            -> ../lib/verkko/scripts/scaffolding/match_graph.py \
                scripts/scaffolding/path_storage.py           -> ../lib/verkko/scripts/scaffolding/path_storage.py \
                scripts/scaffolding/scaffold_graph.py         -> ../lib/verkko/scripts/scaffolding/scaffold_graph.py \
                \
                Snakefile                                     -> ../lib/verkko/Snakefile \
                Snakefiles/1-buildGraph.sm                    -> ../lib/verkko/Snakefiles/1-buildGraph.sm \
                Snakefiles/2-processGraph.sm                  -> ../lib/verkko/Snakefiles/2-processGraph.sm \
                Snakefiles/3-alignONT.sm                      -> ../lib/verkko/Snakefiles/3-alignONT.sm \
                Snakefiles/3-combineONT.sm                    -> ../lib/verkko/Snakefiles/3-combineONT.sm \
                Snakefiles/3-splitONT.sm                      -> ../lib/verkko/Snakefiles/3-splitONT.sm \
                Snakefiles/3-alignTips.sm                     -> ../lib/verkko/Snakefiles/3-alignTips.sm \
                Snakefiles/4-processONT.sm                    -> ../lib/verkko/Snakefiles/4-processONT.sm \
                Snakefiles/5-untip.sm                         -> ../lib/verkko/Snakefiles/5-untip.sm \
                Snakefiles/6-rukki.sm                         -> ../lib/verkko/Snakefiles/6-rukki.sm \
                Snakefiles/6-layoutContigs.sm                 -> ../lib/verkko/Snakefiles/6-layoutContigs.sm \
                Snakefiles/7-buildPackage.sm                  -> ../lib/verkko/Snakefiles/7-buildPackage.sm \
                Snakefiles/7-combineConsensus.sm              -> ../lib/verkko/Snakefiles/7-combineConsensus.sm \
                Snakefiles/7-extractONT.sm                    -> ../lib/verkko/Snakefiles/7-extractONT.sm \
                Snakefiles/7-generateConsensus.sm             -> ../lib/verkko/Snakefiles/7-generateConsensus.sm \
                Snakefiles/8-hicPipeline.sm                   -> ../lib/verkko/Snakefiles/8-hicPipeline.sm \
                Snakefiles/c1-buildStore.sm                   -> ../lib/verkko/Snakefiles/c1-buildStore.sm \
                Snakefiles/c2-findOverlaps.sm                 -> ../lib/verkko/Snakefiles/c2-findOverlaps.sm \
                Snakefiles/c4-findErrors.sm                   -> ../lib/verkko/Snakefiles/c4-findErrors.sm \
                Snakefiles/functions.sm                       -> ../lib/verkko/Snakefiles/functions.sm \
                \
                profiles/config.yaml                          -> ../lib/verkko/profiles/config.yaml \
                profiles/jobscript.sh                         -> ../lib/verkko/profiles/jobscript.sh \
                profiles/slurm-sge-status.sh                  -> ../lib/verkko/profiles/slurm-sge-status.sh \
                profiles/slurm-sge-submit.sh                  -> ../lib/verkko/profiles/slurm-sge-submit.sh


#
#  Rukki doesn't fit into this build system, and is added as a special case.
#
ifeq (${WITHOUT_RUKKI},)
rukki/target/release/rukki: rukki/Cargo.toml \
                            rukki/src/graph.rs \
                            rukki/src/lib.rs \
                            rukki/src/main.rs \
                            rukki/src/pseudo_hap.rs \
                            rukki/src/trio.rs \
                            rukki/src/trio_walk.rs \
                            rukki/src/graph_algos.rs \
                            rukki/src/graph_algos/dfs.rs \
                            rukki/src/graph_algos/scc.rs \
                            rukki/src/graph_algos/superbubble.rs
	cargo build --release --manifest-path rukki/Cargo.toml

doclean: doclean-rukki

.PHONY: doclean-rukki
doclean-rukki:
	rm -rf rukki/target

EXECUTABLES  += rukki/target/release/rukki -> ../lib/verkko/bin/rukki
endif

#
#  MBG doesn't fit into this build system and is added as a special case.
#  The dependency list comes from the MBG makefile.
#
ifeq (${WITHOUT_MBG},)
MBG_DEPS = $(addprefix MBG/src/,\
            $(patsubst %.o,%.cpp,\
                       MBG.o \
                       fastqloader.h CommonUtils.h MBGCommon.h VectorWithDirection.h FastHasher.h SparseEdgeContainer.h HashList.h UnitigGraph.h BluntGraph.h ReadHelper.h HPCConsensus.h ErrorMaskHelper.h CompressedSequence.h ConsensusMaker.h StringIndex.h LittleBigVector.h MostlySparse2DHashmap.h RankBitvector.h TwobitLittleBigVector.h UnitigResolver.h CumulativeVector.h UnitigHelper.h BigVectorSet.h              Serializer.h DumbSelect.h MsatValueVector.h Node.h KmerMatcher.h \
                       fastqloader.o CommonUtils.o MBGCommon.o                       FastHasher.o SparseEdgeContainer.o HashList.o UnitigGraph.o BluntGraph.o              HPCConsensus.o ErrorMaskHelper.o CompressedSequence.o ConsensusMaker.o StringIndex.o                                           RankBitvector.o                         UnitigResolver.o                    UnitigHelper.o BigVectorSet.o ReadHelper.o Serializer.o DumbSelect.o MsatValueVector.o Node.o KmerMatcher.o))

MBG/bin/MBG: $(MBG_DEPS)
	cd MBG ; $(MAKE) all

doclean: doclean-MBG

.PHONY: doclean-MBG
doclean-mbg:
	cd MBG ; $(MAKE) clean

EXECUTABLES  += MBG/bin/MBG -> ../lib/verkko/bin/MBG
endif

#
#  hifioverlapper doesn't fit into this build system and is added as a special case.
#  The dependency list comes from the hifioverlapper makefile.
#
#  To get the two EXECUTABLE rules to work without a race, we need to build
#  each binary independently, and since they both depend on a shared .a, that
#  needs a rule too.
#
ifeq (${WITHOUT_HIFIOVERLAPPER},)
HIFIOVERLAP_DEPS = $(addprefix hifioverlapper/src/,\
                    $(patsubst %.o,%.cpp,\
                               matchchains_index.o matchchains_matchindex.o \
                               MatchIndex.h MinimizerIterator.h TwobitString.h ReadStorage.h UnitigKmerCorrector.h UnitigStorage.h ReadMatchposStorage.h \
                               MatchIndex.o MinimizerIterator.o TwobitString.o ReadStorage.o UnitigKmerCorrector.o UnitigStorage.o ReadMatchposStorage.o))

hifioverlapper/lib/hifioverlapper.a: $(HIFIOVERLAP_DEPS)
	cd hifioverlapper ; $(MAKE) lib/hifioverlapper.a MBG/lib/mbg.a

hifioverlapper/bin/matchchains_index: hifioverlapper/lib/hifioverlapper.a
	cd hifioverlapper ; $(MAKE) bin/matchchains_index

hifioverlapper/bin/matchchains_matchindex: hifioverlapper/lib/hifioverlapper.a
	cd hifioverlapper ; $(MAKE) bin/matchchains_matchindex

doclean: doclean-hifioverlapper

.PHONY: doclean-hifioverlapper
doclean-hifioverlapper:
	cd hifioverlapper ; $(MAKE) clean

EXECUTABLES  += hifioverlapper/bin/matchchains_index      -> ../lib/verkko/bin/matchchains_index \
                hifioverlapper/bin/matchchains_matchindex -> ../lib/verkko/bin/matchchains_matchindex
endif
