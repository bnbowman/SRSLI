add_library (SRSLI_lib
        ReferenceSet.cpp
        SequenceReader.cpp
	Version.cpp
	)

add_executable (srsli
	main.cpp
    )

target_link_libraries(srsli SRSLI_lib ${SEQAN_LIBRARIES})

set(ALL_EXE_TARGETS srsli)

install(TARGETS ${ALL_EXE_TARGETS} RUNTIME DESTINATION bin)
