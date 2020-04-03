start1 = ARGV[0].to_i
start2 = ARGV[1].to_i
allow_mig = ARGV[2]
args = ARGV[2..6].join(" ")
bindir = File.join(File.dirname(__FILE__), "..", "bin")
outdir = File.join(File.dirname(__FILE__), "..", "..", "..", "speciale", "data")

(1..(start1+start2-2)).each do |coals|
	(0..(start1+start2-coals)).each do |pop1|
	    pop2 = start1+start2-coals-pop1
	    system("#{bindir}/run_pure_cutoff #{pop1} #{pop2} #{pop1+pop2-1} #{args} > #{outdir}/next#{pop1}#{pop2}-#{allow_mig}.tsv")
	    system("#{bindir}/run_left_prob #{start1} #{start2} #{pop1} #{pop2} #{args} > #{outdir}/left#{start1}#{start2}-#{allow_mig}-#{pop1}#{pop2}.tsv")
	end
	system("#{bindir}/run_pure_cutoff #{start1} #{start2} #{start1+start2-coals} #{args} > #{outdir}/first#{start1}#{start2}-#{allow_mig}-#{start1+start2-coals}.tsv")
end
