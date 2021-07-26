bos(a,s,u) = 1-(2a)//(2a+s)
bwn(a,s,u) = 1-(2a)//(2a+s+u)
bst(a,s,u) = bwn(a,s,u) - bos(a,s,u)
rst(a,s,u) = bst(a,s,u) // bwn(a,s,u)

bos(2,0,0)
bwn(2,0,0)
bst(2,0,0)
rst(2,0,0)

bos(2,1,0)
bwn(2,1,0)
bst(2,1,0)
rst(2,1,0)

bos(2,0,1)
bwn(2,0,1)
bst(2,0,1)
rst(2,0,1)

bos(2,1,1)
bwn(2,1,1)
bst(2,1,1)
rst(2,1,1)

