(10..80).step(5) do |i|
	puts(i)

	startt = Process.clock_gettime(Process::CLOCK_MONOTONIC)
	system("../bin/run_build_kingman #{i}")
	endt = Process.clock_gettime(Process::CLOCK_MONOTONIC)
	elapsed = endt - startt
	puts(elapsed)
end
