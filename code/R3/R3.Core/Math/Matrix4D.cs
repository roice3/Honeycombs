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

		public Matrix4D( Vector3D[] rows )
		{
			Initialize();
			for( int i=0; i<4; i++ )
				Data[i] = new double[] { rows[i].X, rows[i].Y, rows[i].Z, rows[i].W };
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

		public Vector3D this[int i]
		{
			get
			{
				return new Vector3D( Data[i] );
			}
			set
			{
				Data[i] = new double[] { value.X, value.Y, value.Z, value.W };
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

		/// <summary>
		/// http://www.euclideanspace.com/maths/algebra/matrix/functions/determinant/fourD/index.htm
		/// </summary>
		public double Determinant
		{
			get
			{
				double det =
				Data[0][3] * Data[1][2] * Data[2][1] * Data[3][0] - Data[0][2] * Data[1][3] * Data[2][1] * Data[3][0] - Data[0][3] * Data[1][1] * Data[2][2] * Data[3][0] + Data[0][1] * Data[1][3] * Data[2][2] * Data[3][0] +
				Data[0][2] * Data[1][1] * Data[2][3] * Data[3][0] - Data[0][1] * Data[1][2] * Data[2][3] * Data[3][0] - Data[0][3] * Data[1][2] * Data[2][0] * Data[3][1] + Data[0][2] * Data[1][3] * Data[2][0] * Data[3][1] +
				Data[0][3] * Data[1][0] * Data[2][2] * Data[3][1] - Data[0][0] * Data[1][3] * Data[2][2] * Data[3][1] - Data[0][2] * Data[1][0] * Data[2][3] * Data[3][1] + Data[0][0] * Data[1][2] * Data[2][3] * Data[3][1] +
				Data[0][3] * Data[1][1] * Data[2][0] * Data[3][2] - Data[0][1] * Data[1][3] * Data[2][0] * Data[3][2] - Data[0][3] * Data[1][0] * Data[2][1] * Data[3][2] + Data[0][0] * Data[1][3] * Data[2][1] * Data[3][2] +
				Data[0][1] * Data[1][0] * Data[2][3] * Data[3][2] - Data[0][0] * Data[1][1] * Data[2][3] * Data[3][2] - Data[0][2] * Data[1][1] * Data[2][0] * Data[3][3] + Data[0][1] * Data[1][2] * Data[2][0] * Data[3][3] +
				Data[0][2] * Data[1][0] * Data[2][1] * Data[3][3] - Data[0][0] * Data[1][2] * Data[2][1] * Data[3][3] - Data[0][1] * Data[1][0] * Data[2][2] * Data[3][3] + Data[0][0] * Data[1][1] * Data[2][2] * Data[3][3];
				return det;
			}
		}

		/// <summary>
		/// Gram-Schmidt orthonormalize
		/// </summary>
		public static Matrix4D GramSchmidt( Matrix4D input )
		{
			Matrix4D result = input;
			for( int i=0; i<4; i++ )
			{
				for( int j=0; j<i; j++ )
				{
					// result[j] is already unit length...
					// result[i] -= (result[i] dot result[j])*result[j]
					Vector3D iVec = result[i];
					Vector3D jVec = result[j];
					iVec -= ( iVec.Dot( jVec ) ) * jVec;
					result[i] = iVec;
				}
				result[i].Normalize();
			}

			return result;
		}

		/// <summary>
		/// Gram-Schmidt orthonormalize
		/// </summary>
		public static Matrix4D GramSchmidt( Matrix4D input,
			Func<Vector3D, Vector3D, double> innerProduct, Func<Vector3D, Vector3D> normalize )
		{
			Matrix4D result = input;
			for( int i=0; i<4; i++ )
			{
				for( int j=i+1; j<4; j++ )
				{
					Vector3D iVec = result[i];
					Vector3D jVec = result[j];
					iVec -= innerProduct( iVec, jVec ) * jVec;
					result[i] = iVec;
				}
				result[i] = normalize( result[i] );
			}

			return result;
		}

		/// <summary>
		/// Rotate a vector with this matrix.
		/// </summary>
		public Vector3D RotateVector( Vector3D input )
		{
			Vector3D result = new Vector3D();
			Vector3D copy = new Vector3D( new double[] { input.X, input.Y, input.Z, input.W } );
			for( int i = 0; i < 4; i++ )
			{
				result[i] =
					copy[0] * this[i, 0] +
					copy[1] * this[i, 1] +
					copy[2] * this[i, 2] +
					copy[3] * this[i, 3];
			}
			return result;
		}

		/// <summary>
		/// Returns a matrix which will rotate in a coordinate plane by an angle in radians.
		/// </summary>
		public static Matrix4D MatrixToRotateinCoordinatePlane( double angle, int c1, int c2 )
		{
			Matrix4D result = Matrix4D.Identity();
			result[c1, c1] =  Math.Cos( angle );
			result[c1, c2] = -Math.Sin( angle );
			result[c2, c1] =  Math.Sin( angle );
			result[c2, c2] =  Math.Cos( angle );
			return result;
		}
	}
}
