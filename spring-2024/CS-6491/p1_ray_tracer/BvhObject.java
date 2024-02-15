import java.util.List;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Objects;

// Bounding volume heirarchy acceleration data structure.
class BvhObject extends Object {
    class Node {
        final BBox bbox;
        final Node left, right;
        final int start, end; // Index range for child objects into `BvhObject.objects`, valid only for leaf nodes

        // Constructor for leaf nodes
        Node(BBox bbox, int start, int end) {
            this.bbox = bbox;
            this.start = start;
            this.end = end;
            this.left = this.right = null;
        }

        // Constructor for internal nodes
        Node(BBox bbox, Node left, Node right) {
            this.bbox = bbox;
            this.left = left;
            this.right = right;
            this.start = this.end = -1;
        }

        boolean isLeaf() { return left == null && right == null; }
    }

  // Higher values means bigger bounding boxes with more objects inside them.
  final static int MaxLeafNodeObjectCount = 4;

  final Node root;
  final Surface surface;
  final List<Object> objects;

  BvhObject(List<Object> objects, Surface surface) {
    this.objects = objects;
    this.surface = surface;
    this.root = build(0, objects.size());
  }

  BBox getBBox() { return root.bbox; }
  Hit raycast(Ray ray) { return raycastNode(root, ray); }

  private Node build(int start, int end) {
    final BBox bbox = getBBox(this.objects.subList(start, end));
    if (end - start <= MaxLeafNodeObjectCount) return new Node(bbox, start, end);

    // Sort this range of objects (in place) based on their center along the split axis.
    final int splitAxis = bbox.maxAxis();
    this.objects.subList(start, end).sort(Comparator.comparingDouble(o -> o.getBBox().center().at(splitAxis)));

    final int mid = start + (end - start)/2;
    return new Node(bbox, build(start, mid), build(mid, end));
  }

  private BBox getBBox(List<Object> objects) { return objects.stream().map(Object::getBBox).reduce(BBox.EMPTY, BBox::union); }

  private Hit raycastNode(Node node, Ray ray) {
    if (node == null || !node.bbox.intersects(ray)) return null;

    if (node.isLeaf()) {
      // Uncomment to draw leaf bounding boxes insteaad of the objects inside them.
      //final Intersection isect = node.bbox.intersect(ray);
      //if (isect != null) return new Hit(isect, surface);
      //return null;

      // Check for intersection with all objects in the leaf node.
      return objects.subList(node.start, node.end).stream()
        .map(o -> o.raycast(ray))
        .filter(Objects::nonNull)
        .min(Comparator.comparingDouble(hit -> hit.t))
        .orElse(null);
    }

    // Recurse children.
    final Hit leftHit = raycastNode(node.left, ray), rightHit = raycastNode(node.right, ray);
    if (leftHit != null && rightHit != null) return leftHit.t < rightHit.t ? leftHit : rightHit;

    return leftHit != null ? leftHit : rightHit;
  }
}
