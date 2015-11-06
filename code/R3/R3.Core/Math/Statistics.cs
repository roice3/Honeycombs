namespace R3.Math
{
	using System.Collections.Generic;
	using System.Linq;

	// http://www.codeproject.com/Articles/42492/Using-LINQ-to-Calculate-Basic-Statistics
	// Add more from there as needed.
	public static class Statistics
	{
		public static double Median( this IEnumerable<double> source )
		{
			return ElementAtPercentage( source, 0.5 );
		}

		/// <summary>
		/// This is like Median, but allows you to grab the element an arbitrary percentage along.
		/// percentage should be between 0 and 1.
		/// </summary>
		public static double ElementAtPercentage( this IEnumerable<double> source, double percentage )
		{
			var sortedList = from number in source
							 orderby number
							 select number;

			int count = sortedList.Count();
			int itemIndex = (int)( (double)count * percentage );
			if( count % 2 == 0 ) // Even number of items. 
				return ( sortedList.ElementAt( itemIndex ) +
						sortedList.ElementAt( itemIndex - 1 ) ) / 2;

			// Odd number of items. 
			return sortedList.ElementAt( itemIndex );
		}
	}
}
