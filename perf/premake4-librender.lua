   -- A project defines one build target
   project "LibRender"
      kind "StaticLib"
      language "C++"

      files { render_sources }

      includedirs {
         "../render",
      }

      -- MacOSX. Guess we use gcc.
      configuration { "macosx", "gmake" }

         if _OPTIONS['with-glsl-opt'] then
            defines { 'LSGL_OPTIMIZE_GLSL' }
         end

         -- MPI
         if _OPTIONS['with-mpi'] then
            defines { "LSGL_ENABLE_MPI" }
         end

         -- gcc openmp
         if _OPTIONS['with-openmp'] then
            if _OPTIONS.platform == 'k-cross' or _OPTIONS.platform == 'k-native' then
               buildoptions { "-Kopenmp" }
               linkoptions { "-Kopenmp" }
            else
               buildoptions { "-fopenmp" }
               linkoptions { "-fopenmp" }
            end
         end


      -- Linux specific
      configuration { "linux" }

         if _OPTIONS['with-glsl-opt'] then
            defines { 'LSGL_OPTIMIZE_GLSL' }
            if _OPTIONS.platform == 'k-cross' or _OPTIONS.platform == 'k-native' then
               buildoptions { '-x100', '-Kilfunc', '-Klib' }
            end
         end

         if _OPTIONS['with-k-profile'] then
            defines { "LSGL_ENABLE_K_PROFILE" }
         end

         defines { '__STDC_CONSTANT_MACROS', '__STDC_LIMIT_MACROS' } -- c99
         defines { '_LARGEFILE_SOURCE', '_FILE_OFFSET_BITS=64' }

         -- MPI
         if _OPTIONS['with-mpi'] then
            defines { "LSGL_ENABLE_MPI" }
         end

         if _OPTIONS.platform == 'k-cross' or _OPTIONS.platform == 'k-native' then

            buildoptions { "-Xg -KPIC" }

            -- fj openmp
            if _OPTIONS['with-openmp'] then
               buildoptions { "-Kopenmp" }
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
         if _OPTIONS.platform == 'k-cross' or _OPTIONS.platform == 'k-native' then
            flags { "Symbols" }
         else
            flags { "Symbols" }
         end

         targetname "librender-debug"

      configuration "Release"
         -- defines { "NDEBUG" } -- -NDEBUG

         if _OPTIONS.platform == 'k-cross' or _OPTIONS.platform == 'k-native' then
            -- '-g' invoke -O0 optimization in FCCpx
            buildoptions { "-Kfast" }
         else
            flags { "Symbols", "Optimize", "EnableSSE2" }
         end

         targetname "librender"

