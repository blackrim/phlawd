set (PHLAWD_SOURCES
	DBSeq.cpp
	fasta_util.cpp
	main.cpp
	node.cpp
	sequence.cpp
	SmithWatermanGotoh.cpp
	SQLiteConstructor.cpp
	SQLiteDBController.cpp
	SQLiteProfiler.cpp
	SQLiteTreeNameConvertor.cpp
	SWPS3_debug.cpp
	SWPS3_DynProgr_scalar.cpp
	SWPS3_DynProgr_sse_byte.cpp
	SWPS3_DynProgr_sse_short.cpp
	SWPS3_fasta.cpp
	SWPS3_matrix.cpp
	SWPS3_swps3.cpp
	tree.cpp
	tree_reader.cpp
	tree_utils.cpp
	utils.cpp
)

add_definitions(-msse3 -msse2 -msse)

add_executable (phlawd ${PHLAWD_SOURCES})

target_link_libraries(phlawd ${SQLITE3_LDFLAGS} ${LIBS} ${PROJECT_SOURCE_DIR}/deps/linux/libsqlitewrapped.a)

