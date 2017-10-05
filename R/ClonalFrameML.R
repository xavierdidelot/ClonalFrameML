ClonalFrameML = function(arguments) {
	args = as.character(arguments)
	.C("ClonalFrameMLRmain",PACKAGE="ClonalFrameML",as.integer(length(args)),as.integer(max(nchar(args))),args,EXIT_STATUS=integer(1))$EXIT_STATUS
}
