if (NOT fmt_INCLUDE_DIR)
  set(fmt_INCLUDE_DIR $ENV{fmt_INCLUDE_DIR})
endif ()
if (NOT fmt_LIBRARY)
  set(fmt_LIBRARY $ENV{fmt_LIBRARY})
endif ()

find_path(fmt_INCLUDE_DIR fmt/format.h)
find_library(fmt_LIBRARY NAMES fmt)

mark_as_advanced(fmt_INCLUDE_DIR fmt_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(fmt
    REQUIRED_VARS fmt_INCLUDE_DIR fmt_LIBRARY
)

if(fmt_FOUND AND NOT TARGET fmt::fmt)
    add_library(fmt::fmt UNKNOWN IMPORTED)
    set_target_properties(fmt::fmt PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${fmt_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${fmt_INCLUDE_DIR}"
    )
endif()