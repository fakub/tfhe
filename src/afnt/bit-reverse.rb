def bit_reverse(x, nu)
    r = 0

    nu.times do |i|
        r = (r << 1) | (x & 1)
        x >>= 1
    end

    r
end

nu = 5 ; (1 << (nu-1)).times{|x| printf("%3d | %3d   %s   %3d | %3d\n", 2*x, 2*x+1, (x == (1 << (nu-2)) ? "=>" : "  "), bit_reverse(2*x, nu), bit
_reverse(2*x+1, nu)) }
