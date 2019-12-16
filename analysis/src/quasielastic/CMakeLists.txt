find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/analysis/include ${MYSQL_INCLUDE_DIR} )

include_directories(${INCLUDES})
add_executable(eepler eepler.cpp )
target_link_libraries(eepler ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS eepler
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )

add_executable(incl_ler incl_ler.cpp )
target_link_libraries(incl_ler ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS incl_ler
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )
