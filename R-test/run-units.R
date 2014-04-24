library('RUnit')
library('extraTrees')

test.suite <- defineTestSuite("extraTrees",
                              dirs = file.path("units"),
                              testFileRegexp = '\\.R$')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)

