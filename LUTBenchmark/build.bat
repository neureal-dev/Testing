cmake . -DCMAKE_GENERATOR_PLATFORM=x64 -DCMAKE_TOOLCHAIN_FILE=D:\git\vcpkg\scripts\buildsystems\vcpkg.cmake 
msbuild Project.sln /p:Configuration="Release"