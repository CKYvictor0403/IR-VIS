cmake_minimum_required(VERSION 3.20)

project(FinalProject LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# MSVC: 避免中文訊息編碼警告，統一使用 UTF-8
if (MSVC)
    add_compile_options(/utf-8)
endif()

add_subdirectory(external/wavelib)

# 把 src 底下所有 .cpp 加進來
file(GLOB SRC_FILES
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
)

# 建一個可執行檔 target
add_executable(FinalProject ${SRC_FILES})


# 告訴編譯器 include 在哪裡
target_include_directories(FinalProject
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

# 告訴 FinalProject 去包含 wavelib 的 header
target_include_directories(FinalProject
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/external/wavelib/header
)

# 把 wavelib 這個 target link 進來
target_link_libraries(FinalProject
    PRIVATE
        wavelib
)

# 在 build/ 下建立 output 資料夾
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)