cmake_minimum_required(VERSION 3.20)

project(FinalProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


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

# 在 build/ 下建立 output 資料夾
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)