import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap
import kotlinx.coroutines.experimental.CommonPool
import kotlinx.coroutines.experimental.Deferred
import kotlinx.coroutines.experimental.async
import kotlinx.coroutines.experimental.runBlocking
import java.io.*
import java.nio.charset.StandardCharsets
import java.util.*
import kotlin.experimental.and

object KNucleotide {
	internal val codes = byteArrayOf(-1, 0, -1, 1, 3, -1, -1, 2)
	internal val nucleotides = charArrayOf('A', 'C', 'G', 'T')
	
	internal class Result(var keyLength: Int) {
		var map = Long2IntOpenHashMap()
	}
	
	internal fun sumTwoMaps(map1: Result, map2: Result): Result {
		map2.map.forEach { key, value -> map1.map.addTo(key!!, value!!) }
		return map1
	}
	
	internal fun writeFrequencies(totalCount: Float, frequencies: Result): String {
		val freq = ArrayList<Pair<String, Int>>(frequencies.map.size)
		frequencies.map.forEach { key, cnt ->
			freq.add(Pair(keyToString(key!!, frequencies.keyLength), cnt))
		}
		freq.sortByDescending { it.second }
		val sb = StringBuilder()
		for (entry in freq) {
			sb.append(String.format(Locale.ENGLISH, "%s %.3f\n", entry.first,
															entry.second * 100.0f / totalCount))
		}
		return sb.append('\n').toString()
	}
	
	@Throws(Exception::class)
	internal fun writeCount(futures: List<Deferred<Result>>, nucleotideFragment: String): String = runBlocking(CommonPool) {
		val key = toCodes(nucleotideFragment.toByteArray(StandardCharsets.ISO_8859_1),
											nucleotideFragment.length)
		val k = getKey(key, 0, nucleotideFragment.length)
		var count = 0
		for (future in futures) {
			val f = future.await()
			if (f.keyLength == nucleotideFragment.length) {
				count += f.map.get(k)
			}
		}
		
		count.toString() + "\t" + nucleotideFragment + '\n'
	}
	
	/**
	 * Convert long key to the nucleotides string
	 */
	internal fun keyToString(k: Long, length: Int): String {
		var key = k
		val res = CharArray(length)
		for (i in 0..length - 1) {
			res[length - i - 1] = nucleotides[(key and 0x3).toInt()]
			key = key shr 2
		}
		return String(res)
	}
	
	/**
	 * Get the long key for given byte array of codes at given offset and length
	 * (length must be less than 32)
	 */
	internal fun getKey(arr: ByteArray, offset: Int, length: Int): Long {
		var key: Long = 0
		for (i in offset..offset + length - 1) {
			key = key * 4 + arr[i]
		}
		return key
	}
	
	/**
	 * Convert given byte array (limiting to given length) containing acgtACGT
	 * to codes (0 = A, 1 = C, 2 = G, 3 = T) and returns new array
	 */
	internal fun toCodes(sequence: ByteArray, length: Int): ByteArray {
		val result = ByteArray(length)
		for (i in 0..length - 1) {
			result[i] = codes[(sequence[i] and 0x7).toInt()]
		}
		return result
	}
	
	@Throws(IOException::class)
	fun read(`is`: InputStream): ByteArray {
		var line: String
		val `in` = BufferedReader(InputStreamReader(`is`, StandardCharsets.ISO_8859_1))
		line = `in`.readLine()
		while (line != null) {
			if (line.startsWith(">THREE"))
				break
			line = `in`.readLine()
		}
		
		var bytes = ByteArray(1048576)
		var position = 0
		line = `in`.readLine()
		while (line.isNotEmpty() && line[0] != '>') {
			if (line.length + position > bytes.size) {
				val newBytes = ByteArray(bytes.size * 2)
				System.arraycopy(bytes, 0, newBytes, 0, position)
				bytes = newBytes
			}
			for (i in 0..line.length - 1)
				bytes[position++] = line[i].toByte()
			line = `in`.readLine() ?: ""
		}
		
		return toCodes(bytes, position)
	}
	
	internal fun createFragmentMap(sequence: ByteArray, offset: Int, fragmentLength: Int): Deferred<Result> = async(CommonPool) {
		val res = Result(fragmentLength)
		val map = res.map
		val lastIndex = sequence.size - fragmentLength + 1
		var index = offset
		while (index < lastIndex) {
			map.addTo(getKey(sequence, index, fragmentLength), 1)
			index += fragmentLength
		}
		
		res
	}
	
	internal fun createFragmentTasks(sequence: ByteArray,
																	 fragmentLengths: IntArray): ArrayList<Deferred<Result>> {
		val tasks = ArrayList<Deferred<Result>>()
		for (fragmentLength in fragmentLengths) {
			(0..fragmentLength - 1)
				.mapTo(tasks) {  createFragmentMap(sequence, it, fragmentLength) }
		}
		return tasks
	}
}

fun main(args: Array<String>) = runBlocking(CommonPool) {
	val sequence = KNucleotide.read(
		FileInputStream("${System.getProperty("user.dir")}/src/main/resources/input_25kk.txt"))
	val begin = System.nanoTime()
	
	val fragmentLengths = intArrayOf(1, 2, 3, 4, 6, 12, 18)
	val futures = KNucleotide.createFragmentTasks(sequence, fragmentLengths)
	
	val sb = StringBuilder()
	
	sb.append(KNucleotide.writeFrequencies(sequence.size.toFloat(), futures[0].await()))
	sb.append(KNucleotide.writeFrequencies((sequence.size - 1).toFloat(),
																				 KNucleotide.sumTwoMaps(futures[1].await(), futures[2].await())))
	
	val nucleotideFragments = arrayOf("GGT", "GGTA", "GGTATT", "GGTATTTTAATT", "GGTATTTTAATTTATAGT")
	for (nucleotideFragment in nucleotideFragments) {
		sb.append(KNucleotide.writeCount(futures, nucleotideFragment))
	}
	
	print(sb)
	println((System.nanoTime() - begin).toDouble() / 1000_000_000.toDouble())
}

