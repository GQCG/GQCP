import os
import glob


def prepend_code(text, prepending_characters):
    """Prepend and :return :param text (a string of source code) with :param prepending_characters
    at the beginning of a new line.

    This function is used to insert commenting characters in front of a block of text.
    """
    text = prepending_characters + text  # prepend the first line
    text = text.replace('\n', '\n' + prepending_characters)  # prepend the other lines
    return text


# Find the relevant files and directories
tools_directory_path = os.path.dirname(os.path.realpath(__file__))
header_path = os.path.join(tools_directory_path, 'HEADER')
library_directory_path = os.path.dirname(tools_directory_path)  # parent of 'tools'

include_directory_path = os.path.join(library_directory_path, 'include')
header_files = glob.glob(include_directory_path + "/*.hpp")

src_directory_path = os.path.join(library_directory_path, 'src')
source_files = glob.glob(src_directory_path + "/*.cpp")

exe_directory_path = os.path.join(library_directory_path, 'exe')
executable_files = glob.glob(exe_directory_path + "/*.cpp")


# Add or update the header in every header (.hpp) or source file (.cpp)
with open(header_path, 'r') as f_header:
    header_text = f_header.read()  # read in the whole header
    header_text_commented = prepend_code(header_text, '// ')  # '//' for comments in C++
    header_text_commented += '\n'  # we need an extra (non-commented) newline at the end

    for filename in (header_files + source_files + executable_files):

        with open(filename, 'r') as f_original:
            original_source_code = f_original.read()  # read the whole source code

        # Add the header in the case that the copyright header isn't included yet
        if original_source_code[0] == '#':
            file_is_raw = True
        elif original_source_code[0] == '/':
            file_is_raw = False
        else:
            raise ValueError("{filename} has an unexpected character {character} at position 0".
                             format(filename=filename, character=original_source_code[0]))

        with open(filename, 'w') as f_modified:
            if file_is_raw:  # needs 'adding' of the copyright header
                    f_modified.write(header_text_commented + original_source_code)

            else:  # needs 'updating' of the copyright header
                # Splitting at the first occurrence splits off the copyright header, but we'll have
                # to add the '#' again later
                original_source_code_without_copyright = original_source_code.split('#', 1)[1]

                f_modified.write(header_text_commented + '#'
                                 + original_source_code_without_copyright)
