cmake_minimum_required(VERSION 3.15)
project(star_apparent_catalog LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ERFA 子集（原样收录, LGPL-3.0, 见 third_party/erfa/LICENSE）
set(ERFA_SRC
    third_party/erfa/epv00.c
    third_party/erfa/pmpx.c
    third_party/erfa/ab.c
    third_party/erfa/ld.c
    third_party/erfa/pdp.c
    third_party/erfa/pxp.c
    third_party/erfa/pn.c
    third_party/erfa/anp.c
    third_party/erfa/dat.c
    third_party/erfa/erfadatextra.c
    third_party/erfa/zp.c
    third_party/erfa/pm.c
    third_party/erfa/sxp.c
    third_party/erfa/cal2jd.c
    third_party/erfa/jd2cal.c
)

add_library(apparent_catalog STATIC
    apparent_catalog.c
    ${ERFA_SRC}
)
target_include_directories(apparent_catalog
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/third_party/erfa
)

if(MSVC)
    # 源码为 UTF-8（无 BOM）编码, 必须显式声明
    target_compile_options(apparent_catalog PRIVATE /utf-8)
    target_compile_options(apparent_catalog PRIVATE /W4)
    target_compile_definitions(apparent_catalog PRIVATE _CRT_SECURE_NO_WARNINGS)
    # ERFA 第三方源码不做高级别告警检查
    set_source_files_properties(${ERFA_SRC} PROPERTIES COMPILE_OPTIONS "/W3")
else()
    target_compile_options(apparent_catalog PRIVATE -Wall -Wextra)
    if(UNIX)
        target_link_libraries(apparent_catalog PUBLIC m)
    endif()
endif()

add_executable(test_apparent_catalog test_apparent_catalog.c)
target_link_libraries(test_apparent_catalog PRIVATE apparent_catalog)
if(MSVC)
    target_compile_options(test_apparent_catalog PRIVATE /utf-8 /W4)
    target_compile_definitions(test_apparent_catalog PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()
