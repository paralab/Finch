this_dir=$(pwd)
exec_dir="$(dirname "$this_dir")"
mpirun_command=$(which mpirun)

echo $this_dir
echo $exec_dir
echo $mpirun_command



# run every regression test
echo ---------------------- Test ---------------------
cd $this_dir
if [ $# -eq 0 ]
  then
    # compile everything (assume already cmaked)
    echo ---------------------- Compile ---------------------
    # cd $exec_dir/$build_type/examples/Basic
    for build_folder in $exec_dir/cmake-build*
    do
      cd $build_folder
      make -j 4 &
      cd $this_dir
    done
    wait
    echo "build all complete"
    for sub_dir in */ ;
    do
      cd $sub_dir
      echo ---------------------- Testing $sub_dir ---------------------
      python3 *_regression*.py --runner=$mpirun_command
      echo -------------------------------------------------------------
      cd $this_dir
    done
  else
    echo ---------------------- Compile ---------------------
    # cd $exec_dir/$build_type/examples/Basic
    for build_folder in $exec_dir/cmake-build*/examples/${1}
    do
      cd $build_folder
      make -j 4 &
      cd $this_dir
    done
    wait
    for sub_dir in $this_dir/${1}/ ;
      do
        cd $sub_dir
        echo ---------------------- Testing $sub_dir ---------------------
        python3 *_regression*.py --runner=$mpirun_command
        echo -------------------------------------------------------------
        cd $this_dir
      done
fi
