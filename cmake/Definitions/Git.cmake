# Get git-related info


# Get the long git SHA1 (https://stackoverflow.com/a/21028226/7930415)
execute_process(COMMAND "${GIT_EXECUTABLE}" describe --always --abbrev=40 --dirty
                WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                OUTPUT_VARIABLE GIT_SHA1
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)
