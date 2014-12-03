module TestQuadrotor

importall BatchTest
importall StandardQuadrotor

export run_tests, load_test, series, file_name
export data_series, test_types

# file_name
# system_type dim dim dim.test
data_series = {
  "1d" => [1, 1, 2].^([0:10]'),
  "2d" => [1, 2, 2].^([0:6]'),
  "3d" => [2, 2, 2].^([0:4]')
  }

test_types = [multi_agent_system, relative_multi_agent_system, payload_system,
  relative_payload_system]

function file_name(system_fun, series, index::Integer)
  dim = data_series[series][:,index]
  file_name(system_fun, series, dim)
end

function file_name(system_fun, series, dim::Vector)
  x = dim[1]
  y = dim[2]
  z = dim[3]
  "./test/$system_fun#$series#$x\_$y\_$z"
end

function create_test(system_fun, series, index)
  dim = data_series[series][:, index]
  Test(file_name(system_fun, series, dim),
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

function series(system_fun, series_name)
  tests = data_series[series_name]
  [(tests[:,1], file_name(system_fun, series_name, tests[:,1])) for i = 1:size(tests,2)]
end

function load_test(fun::Integer, series::Integer, index::Integer)
  series_names = collect(keys(data_series))
  series_name = series_names[series]
  name = file_name(test_types[fun], series_name, index)
  load_test(name)
end

end
