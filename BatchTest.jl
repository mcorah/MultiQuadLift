module BatchTest

using HDF5, JLD

export Test, execute_test

type Test
  file_name::String
  result_function::Function
end

function execute_test(test::Test; overwrite = false)
  file_name = test.file_name
  if overwrite == true || filemode(file_name) == 0x0
    println("Running $file_name")
    @time result = test.result_function()
    println("Saving $file_name")
    @save file_name result
  else
    println("Skipping: $file_name already exists")
  end
  println()
end

load_test(test::Test) = load_test(test.file_name) 
function load_test(file_name)
  @load file_name result
end

end
