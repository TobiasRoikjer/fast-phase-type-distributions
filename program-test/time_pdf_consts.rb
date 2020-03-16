(10..60).step(5) do |i|
	puts(i)

	res = `../bin/run_time pdf_constants #{i} 1`
	puts(res.split(":")[1])
end
