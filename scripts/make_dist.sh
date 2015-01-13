#!/bin/sh
echo "Make dist dir"
mkdir dist
echo "Copy gcc libraries"
cp /usr/local/lib/libstdc++.6.dylib dist/
#cp /usr/local/lib/libgomp.dylib dist/
cp /usr/local/lib/libgcc_s.1.dylib dist/
echo "Change permission of libgcc_s"
chmod +wx dist/libgcc_s.1.dylib

echo "Copy krender/KAnim"
cp app/KAnimRender/krender dist/
cp app/KAnimRender/KAnim dist/

echo "Copy libLSGLES"
cp gles/libLSGLES.dylib dist/

echo "Change id of gcc libs."
cd dist
install_name_tool -id ./libstdc++.6.dylib ./libstdc++.6.dylib
install_name_tool -id ./libgcc_s.1.dylib ./libgcc_s.1.dylib
#install_name_tool -id ./libgomp.dylib ./libgomp.dylib

echo "Resolve dependencies in libstd/libgomp"
install_name_tool -change /usr/local/lib/libgcc_s.1.dylib ./libgcc_s.1.dylib libstdc++.6.dylib
#install_name_tool -change /usr/local/lib/libgcc_s.1.dylib ./libgcc_s.1.dylib libgomp.dylib

echo "Resolve dependencies in libLSGLES.dylib"
install_name_tool -change /usr/local/lib/libstdc++.6.dylib ./libstdc++.6.dylib libLSGLES.dylib
install_name_tool -change /usr/local/lib/libgcc_s.1.dylib ./libgcc_s.1.dylib libLSGLES.dylib

echo "Copy shaders"
cd ../
cp -Rf app/KAnimRender/krender_node/shader dist/

install_name_tool -change ./libLSGLES.dylib @executable_path/libLSGLES.dylib dist/KAnim


echo "Done"
