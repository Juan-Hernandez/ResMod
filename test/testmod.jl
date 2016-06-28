module testmod
	export detsvector, testX, testfun!
	testX=Array(Float64,3,3,40,40)
	testX[:]=1:14400
	testX[1,1,:,:]=-1.0
	detsvector=Array(Float64,1600,1)	
	@time detsvector[:]=pmap(det, [testX[:,:,i] for i=1:1600]) 
""" The previous one is way better, less memory, less time. ONLY sub() those arrays preallocating outputs """
	#@time detsvector=pmap(det, [sub(testX, :, :, i) for i=1:1600])
	#@time map!(det, sub(detsvector, :, 1), [sub(testX, :, :, i) for i=1:16])
	function testfun!(r::Any, x::Int64, y::Int64)
		r[:]=x+y
		nothing
	end
end

