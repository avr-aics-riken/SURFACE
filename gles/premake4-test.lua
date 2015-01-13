   -- A project building unit test
   project "LSGLESTest"
      kind "ConsoleApp"
      language "C++"

      files { test_sources, gtest_sources, gles_sources, render_sources }

      includedirs {
         "./",
         "../render",
         "../glsl",
         "../include",
         "../third_party/gtest-1.7.0/include",
         "../third_party/gtest-1.7.0/"
      }

      if _OPTIONS['with-debugtrace'] then
         defines { "LSGL_DEBUG_TRACE" }
      end

      -- for gtest
      if k_platform then 
         defines { 'GTEST_HAS_TR1_TUPLE=0' }
      end

      -- MacOSX. Guess we use gcc.
      configuration { "macosx", "gmake" }

         -- Screen parallel
         if _OPTIONS['with-screen-parallel'] then
            defines { "LSGL_ENABLE_SCREEN_PARALLEL" }
         end

         -- MPI
         if _OPTIONS['with-mpi'] then
            defines { "LSGL_ENABLE_MPI" }
         end

         -- gcc openmp
         if _OPTIONS['with-openmp'] then
            buildoptions { "-fopenmp" }
            linkoptions { "-fopenmp" }

            -- statically link libgomp.
            linkoptions { "`${CXX} -print-file-name=libgomp.a`" }

            -- statically link libgcc and libstdc++
            linkoptions { "-static-libgcc", "-static-libstdc++" }
         end

      -- Windows specific. Assumes visual studio.
      configuration { "windows" }

         defines { "LSGL_EXPORT" }

         includedirs { './compat' } -- stdint.h 
         defines { 'NOMINMAX', '_LARGEFILE_SOURCE', '_FILE_OFFSET_BITS=64' }
         defines { '__STDC_CONSTANT_MACROS', '__STDC_LIMIT_MACROS' } -- c99

         -- Screen parallel
         if _OPTIONS['with-screen-parallel'] then
            defines { "LSGL_ENABLE_SCREEN_PARALLEL" }
         end

         -- MPI
         if _OPTIONS['with-mpi'] then
            defines { "LSGL_ENABLE_MPI" }
         end

      -- Linux specific
      configuration { "linux" }

         if not k_platform then 
            links { "dl" }
         end

         links { "pthread" }

         if _OPTIONS['with-k-profile'] then
            defines { "LSGL_ENABLE_K_PROFILE" }
         end

         defines { '__STDC_CONSTANT_MACROS', '__STDC_LIMIT_MACROS' } -- c99
         defines { '_LARGEFILE_SOURCE', '_FILE_OFFSET_BITS=64' }

         -- Screen parallel
         if _OPTIONS['with-screen-parallel'] then
            defines { "LSGL_ENABLE_SCREEN_PARALLEL" }
         end

         -- MPI
         if _OPTIONS['with-mpi'] then
            defines { "LSGL_ENABLE_MPI" }
         end

         if k_platform then 
            buildoptions { "-Xg -KPIC" }

            -- fj openmp
            if _OPTIONS['with-openmp'] then
               buildoptions { "-Kopenmp" }
               linkoptions { "-Kopenmp" }
            end
         else
            -- gcc openmp
            if _OPTIONS['with-openmp'] then
               buildoptions { "-fopenmp" }
               linkoptions { "-fopenmp" }
            end
         end

      configuration "Debug"
         defines { "DEBUG" } -- -DDEBUG
         if k_plaform then 
            flags { "Symbols" }
         else
            flags { "Symbols", "EnableSSE2" }
         end

         if _OPTIONS['with-mpi'] then
            if _OPTIONS['with-screen-parallel'] then
               targetname "test_LSGLESd_mpi_screenparallel"
            else
               targetname "test_LSGLESd_mpi"
            end
         else
            targetname "test_LSGLESd"
         end

      configuration "Release"
         if k_platform then 
            -- '-g' invoke -O0 optimization in FCCpx
            buildoptions { "-Kfast" }
         else
            flags { "Symbols", "Optimize", "EnableSSE2" }
         end

         if _OPTIONS['with-mpi'] then
            if _OPTIONS['with-screen-parallel'] then
               targetname "test_LSGLES_mpi_screenparallel"
            else
               targetname "test_LSGLES_mpi"
            end
         else
            targetname "test_LSGLES"
         end
