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

benchmark_sources = {
   "bm-point-bvh.cc",
   }

harness_sources = {
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

function newplatform(plf)
    local name = plf.name
    local description = plf.description

    -- Register new platform
    premake.platforms[name] = {
        cfgsuffix = "_"..name,
        iscrosscompiler = true
    }

    -- Allow use of new platform in --platfroms
    table.insert(premake.option.list["platform"].allowed, { name, description })
    table.insert(premake.fields.platforms.allowed, name)

    -- Add compiler support
    -- gcc
    premake.gcc.platforms[name] = plf.gcc
    --other compilers (?)
end

newplatform {
    name = "k-cross",
    description = "K/FX10 cross compiler",
    gcc = {
        cc = "fccpx",
        cxx = "FCCpx",
        cppflags = "-MMD -Xg"
    }
}

newplatform {
    name = "k-native",
    description = "K/FX10 native compiler",
    gcc = {
        cc = "fcc",
        cxx = "FCC",
        cppflags = "-MMD -Xg"
    }
}

newplatform {
    name = "icc",
    description = "Intel C/C++ compiler",
    gcc = {
        cc = "icc",
        cxx = "icpc",
        cppflags = "-MMD"
    }
}

-- premake4.lua
solution "RenderBenchmarkSolution"
   configurations { "Release", "Debug" }

   -- default to x64 platform if none has been specified
   _OPTIONS.platform = _OPTIONS.platform or 'x64'

   -- overwrite the native platform with the options::platform
   premake.gcc.platforms['Native'] = premake.gcc.platforms[_OPTIONS.platform]

   -- A project defines one build target
   project "RenderBenchmark"
      kind "ConsoleApp"
      language "C++"

      -- dependency
      links "LibRender"

      files { libbenchmark_sources, benchmark_sources, harness_sources }

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

         if _OPTIONS['with-openmp'] then
            if (_OPTIONS.platform == 'k-cross' or _OPTIONS.platform == 'k-native') then
               buildoptions { "-Kopenmp" }
               linkoptions { "-Kopenmp" }
            else
               buildoptions { "-fopenmp" }
               linkoptions { "-fopenmp" }
            end
         end

      -- Linux specific
      configuration { "linux" }

         links { "pthread", "rt" }

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

         if _OPTIONS['with-openmp'] then
            if (_OPTIONS.platform == 'k-cross' or _OPTIONS.platform == 'k-native') then
               buildoptions { "-Kopenmp" }
               linkoptions { "-Kopenmp" }
            else
               buildoptions { "-fopenmp" }
               linkoptions { "-fopenmp" }
            end
         end

      configuration "Debug"
         defines { "DEBUG" } -- -DDEBUG

         targetname "render-benchmark-debug"

      configuration "Release"

         flags { "Symbols", "Optimize" }

         targetname "render-benchmark"

      -- Include lib
      dofile "premake4-librender.lua"
