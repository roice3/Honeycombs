namespace R3.Math
{
	using R3.Geometry;
	using System;
	using Math = System.Math;

	public class Matrix4D
	{
		public Matrix4D()
		{
			Initialize();
		}

		public Matrix4D( double[,] data )
		{
			Initialize();
			for( int i=0; i<4; i++ )
			for( int j=0; j<4; j++ )
				Data[i][j] = data[i,j];
		}

		private void Initialize()
		{
			Data = new double[4][];
			for( int i=0; i<4; i++ )
				Data[i] = new double[4];
		}

		public double[][] Data { get; set; }

		public Matrix4D Clone()
		{
			Matrix4D result = new Matrix4D();
			for( int i=0; i<4; i++ )
			for( int j=0; j<4; j++ )
				result.Data[i][j] = this.Data[i][j];
			return result;
		}

		public static Matrix4D Identity()
		{
			Matrix4D result = new Matrix4D();
			for( int i=0; i<4; i++ )
				result[i,i] = 1;
			return result;
		}

		/// <summary>
		/// Mixing multidim and jagged array notation here, but whatevs.
		/// </summary>
		public double this[int i, int j]
		{
			get
			{
				return Data[i][j];
			}
			set
			{
				Data[i][j] = value;
			}
		}

		public VectorND this[int i]
		{
			get
			{
				return new VectorND( Data[i] );
			}
			set
			{
				Data[i] = value.X;
			}
		}

		public static Matrix4D operator +( Matrix4D m1, Matrix4D m2 )
		{
			Matrix4D result = new Matrix4D();
			for( int i=0; i<4; i++ )
			for( int j=0; j<4; j++ )
				result[i, j] = m1[i, j] + m2[i, j];
			return result;
		}

		public static Matrix4D operator *( Matrix4D m1, Matrix4D m2 )
		{
			Matrix4D result = new Matrix4D();
			for( int i=0; i<4; i++ )
			for( int j=0; j<4; j++ )
			for( int k=0; k<4; k++ )
				result[i, j] += m1[i, k] * m2[k, j];
			return result;
		}

		public static Matrix4D operator *( Matrix4D m, double s )
		{
			Matrix4D result = new Matrix4D();
			for( int i=0; i<4; i++ )
			for( int j=0; j<4; j++ )
				result[i, j] = m[i, j] * s;
			return result;
		}

		public static Matrix4D Transpose( Matrix4D m )
		{
			Matrix4D result = new Matrix4D();
			for( int i=0; i<4; i++ )
			for( int j=0; j<4; j++ )
				result[i, j] = m[j, i];
			return result;
		}
	}
}
