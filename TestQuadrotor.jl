module TestQuadrotor

using BatchTest
using StandardQuadrotor

export run_tests

# file_name
# system_type dim dim dim.test
data_series = {
  "1d_z" => [1, 1, 2].^([0:10]'),
  "1d_y" => [1, 2, 1].^([0:10]'),
  "2d_z" => [2, 2, 1].^([0:8]'),
  "2d_y" => [2, 1, 2].^([0:8]'),
  "3d" => [2, 2, 2].^([0:5]')
  }

test_types = [multi_agent_system, relative_multi_agent_system, payload_system,
  relative_payload_system]

file_name(system_fun, series, x, y, z) = "./test/$system_fun#$series#$x\_$y\_$z"

function create_test(system_fun, series, index)
  dim = data_series[series][:, index]
  x = dim[1]
  y = dim[2]
  z = dim[3]
  Test(file_name(system_fun, series, x, y, z),
    ()->run_test(system_fun, dim))
end

tests = vcat([ 
    vcat([ 
      [ create_test(system_fun, series, index) for system_fun=test_types]
    for index = 1:size(data_series[series],2)]...)
  for series = keys(data_series)]...)

function run_tests(overwrite = false)
  println("Running Tests\n")
  for test = tests
    execute_test(test, overwrite = overwrite)
  end
end

end
