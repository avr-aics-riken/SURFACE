-- OpenMP
newoption {
   trigger     = "with-openmp",
   description = "Use OpenMP."
}

-- MPI
newoption {
   trigger     = "with-mpi",
   description = "Use MPI."
}

-- Tentative
newoption {
   trigger     = "with-glsl-opt",
   description = "Use Optimized GLSL codepath."
}

main_sources = {
   "main.cc",
   "harness-accel-point.cc",
   "harness-random-points.cc",
   "harness-random-triangles.cc",
   "harness-prefix-tree.cc",
   "harness-json.cc",
   }

render_sources = {
   "../render/render_accel_triangle.cc",
   "../render/render_accel_line.cc",
   "../render/render_accel_particle.cc",
   "../render/render_camera.cc",
   "../render/render_prim_triangle.cc",
   "../render/render_texture.cc",
   "../render/render_bvh_tree.cc",
   "../render/render_prefix_tree_util.cc",
   "../render/tinymt64.cpp",
   }

libbenchmark_sources = {
   "./deps/hayai/src/hayai_posix_main.cpp",
}

-- premake4.lua
solution "RenderBenchmarkSolution"
   configurations { "Release", "Debug" }

   if (os.is("windows")) then
      platforms { "x32", "x64", "native" }
   else
      platforms { "native", "x32", "x64" }
   end

   -- A project defines one build target
   project "RenderBenchmark"
      kind "ConsoleApp"
      language "C++"

      -- dependency
      links "LibRender"

      files { libbenchmark_sources, main_sources }

      includedirs {
         "./deps/hayai/src/",
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
            buildoptions { "-fopenmp" }
            linkoptions { "-fopenmp" }
         end

      -- Linux specific
      configuration { "linux" }

         links { "pthread" }

         if _OPTIONS['with-glsl-opt'] then
            defines { 'LSGL_OPTIMIZE_GLSL' }
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

         -- gcc openmp
         if _OPTIONS['with-openmp'] then
            buildoptions { "-fopenmp" }
            linkoptions { "-fopenmp" }
         end

      configuration "Debug"
         defines { "DEBUG" } -- -DDEBUG

         targetname "render-benchmark-debug"

      configuration "Release"

         flags { "Symbols", "Optimize" }

         targetname "render-benchmark"

      -- Include lib
      dofile "premake4-librender.lua"
