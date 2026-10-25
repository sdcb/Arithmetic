<Query Kind="Program">
  <NuGetReference Prerelease="true">Sdcb.Arithmetic</NuGetReference>
  <Namespace>Sdcb.Arithmetic</Namespace>
  <Namespace>System.Diagnostics</Namespace>
  <Namespace>System.Security.Cryptography</Namespace>
</Query>

// Textbook RSA built from BigInt primitives: key generation (NextPrime, ModInverse),
// encryption/decryption (ModPow, with the CRT speed-up) and PKCS#1 v1.5 signatures that
// interoperate with System.Security.Cryptography.RSA in both directions.
// FOR LEARNING ONLY: no constant-time arithmetic, no padding for encryption.
void Main()
{
	int bits = 2048;
	Stopwatch sw = Stopwatch.StartNew();
	RsaKey key = RsaKey.Generate(bits);
	Console.WriteLine($"generated {bits}-bit key in {sw.ElapsedMilliseconds} ms");
	Console.WriteLine($"n = {key.N.ToString(16)[..32]}... ({key.N.GetBitLength()} bits)");

	// --- Encrypt / decrypt -----------------------------------------------------------
	string message = "Hello, Sdcb.Arithmetic!";
	BigInt m = BigInt.FromMagnitudeBytes(Encoding.UTF8.GetBytes(message), bigEndian: true);
	BigInt c = BigInt.ModPow(m, key.E, key.N);
	BigInt plain = BigInt.ModPow(c, key.D, key.N), crt = key.DecryptCrt(c); // warm up
	sw.Restart();
	for (int i = 0; i < 20; i++) plain = BigInt.ModPow(c, key.D, key.N);
	double plainMs = sw.Elapsed.TotalMilliseconds / 20;
	sw.Restart();
	for (int i = 0; i < 20; i++) crt = key.DecryptCrt(c);
	double crtMs = sw.Elapsed.TotalMilliseconds / 20;
	Console.WriteLine($"decrypted: \"{Encoding.UTF8.GetString(plain.ToMagnitudeBytes(bigEndian: true))}\"");
	Console.WriteLine($"decrypt: c^d mod n {plainMs:F2} ms, CRT {crtMs:F2} ms (same result: {crt == plain})");

	// --- Signatures, cross-checked with .NET's RSA -----------------------------------
	using RSA rsa = RSA.Create();
	rsa.ImportParameters(key.ToRsaParameters());
	byte[] data = Encoding.UTF8.GetBytes("sign me");

	byte[] ours = key.Sign(data);
	byte[] theirs = rsa.SignData(data, HashAlgorithmName.SHA256, RSASignaturePadding.Pkcs1);
	Console.WriteLine($".NET verifies our signature: {rsa.VerifyData(data, ours, HashAlgorithmName.SHA256, RSASignaturePadding.Pkcs1)}");
	Console.WriteLine($"we verify .NET's signature : {key.Verify(data, theirs)}");
	Console.WriteLine($"signatures identical       : {ours.AsSpan().SequenceEqual(theirs)}");
}

record RsaKey(int Bits, BigInt N, BigInt E, BigInt D, BigInt P, BigInt Q, BigInt DP, BigInt DQ, BigInt QInv)
{
	int Length => (Bits + 7) / 8;

	public static RsaKey Generate(int bits)
	{
		BigInt e = 65537;
		BigInt p, q, n;
		do
		{
			p = RandomPrime(bits / 2);
			q = RandomPrime(bits - bits / 2);
			n = p * q;
		} while (p == q || n.GetBitLength() != bits || BigInt.Gcd(e, (p - 1) * (q - 1)) != 1);
		BigInt d = BigInt.ModInverse(e, BigInt.Lcm(p - 1, q - 1));
		return new RsaKey(bits, n, e, d, p, q, d % (p - 1), d % (q - 1), BigInt.ModInverse(q, p));
	}

	// A random prime with exactly `bits` bits (top two bits set so p*q has the full length).
	static BigInt RandomPrime(int bits)
	{
		byte[] buf = RandomNumberGenerator.GetBytes((bits + 7) / 8);
		BigInt x = BigInt.FromMagnitudeBytes(buf) >> (buf.Length * 8 - bits);
		x |= BigInt.One << (bits - 1) | BigInt.One << (bits - 2);
		return x.NextPrime();
	}

	public BigInt DecryptCrt(BigInt c)
	{
		BigInt m1 = BigInt.ModPow(c, DP, P), m2 = BigInt.ModPow(c, DQ, Q);
		BigInt h = QInv * (m1 - m2) % P;
		if (h.Sign < 0) h += P;
		return m2 + h * Q;
	}

	// EMSA-PKCS1-v1_5 with SHA-256: 00 01 FF..FF 00 || DigestInfo || hash.
	BigInt EncodePkcs1(byte[] data)
	{
		byte[] prefix = Convert.FromHexString("3031300d060960864801650304020105000420");
		byte[] em = new byte[Length];
		byte[] hash = SHA256.HashData(data);
		em[1] = 0x01;
		int tail = prefix.Length + hash.Length;
		em.AsSpan(2, em.Length - 3 - tail).Fill(0xFF);
		prefix.CopyTo(em, em.Length - tail);
		hash.CopyTo(em, em.Length - hash.Length);
		return BigInt.FromMagnitudeBytes(em, bigEndian: true);
	}

	public byte[] Sign(byte[] data) => Fixed(DecryptCrt(EncodePkcs1(data)), Length);

	public bool Verify(byte[] data, byte[] signature) =>
		BigInt.ModPow(BigInt.FromMagnitudeBytes(signature, bigEndian: true), E, N) == EncodePkcs1(data);

	public RSAParameters ToRsaParameters()
	{
		int half = (Length + 1) / 2;
		return new RSAParameters
		{
			Modulus = Fixed(N, Length), Exponent = E.ToMagnitudeBytes(bigEndian: true), D = Fixed(D, Length),
			P = Fixed(P, half), Q = Fixed(Q, half), DP = Fixed(DP, half), DQ = Fixed(DQ, half), InverseQ = Fixed(QInv, half),
		};
	}

	// Big-endian magnitude left-padded to a fixed length.
	static byte[] Fixed(BigInt v, int length)
	{
		byte[] raw = v.ToMagnitudeBytes(bigEndian: true);
		byte[] result = new byte[length];
		raw.CopyTo(result, length - raw.Length);
		return result;
	}
}
