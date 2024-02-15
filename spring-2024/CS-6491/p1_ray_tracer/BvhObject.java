import java.util.List;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Objects;

// Bounding volume heirarchy acceleration data structure.
class BvhObject extends Object {
  class Node {
    final BBox bbox;
    final Node left, right;
    final List<Object> objects; // Non-null only for leaf nodes

    // For leaf node
    Node(BBox bbox, List<Object> objects) {
      this.bbox = bbox;
      this.objects = objects;
      this.left = null;
      this.right = null;
    }

    // For internal node
    Node(BBox bbox, Node left, Node right) {
      this.bbox = bbox;
      this.left = left;
      this.right = right;
      this.objects = null;
    }

    boolean isLeaf() { return objects != null; }
  }

  final static int MaxLeafNodeObjectCount = 6;

  final Node root;
  final Surface surface;

  BvhObject(List<Object> objects, Surface surface) {
    this.surface = surface;
    this.root = build(objects, 0, objects.size());
  }

  BBox getBBox() { return root.bbox; }
  Hit raycast(Ray ray) { return raycastNode(root, ray); }

  private Node build(List<Object> objects, int start, int end) {
    if (end - start <= MaxLeafNodeObjectCount) {
      final List<Object> leafObjects = new ArrayList(objects.subList(start, end));
      return new Node(getBBox(leafObjects), leafObjects);
    }

    // Simple split criterion: median of the bounding box's longest axis
    final BBox bbox = getBBox(objects);
    final int splitAxis = bbox.maxAxis();
    objects.sort(Comparator.comparingDouble(o -> o.getBBox().center().at(splitAxis)));
    final int mid = start + (end - start) / 2;
    return new Node(bbox, build(objects, start, mid), build(objects, mid, end));
  }

  private BBox getBBox(List<Object> objects) {
    final float MAX = Float.MAX_VALUE;
    return objects.stream()
      .map(Object::getBBox)
      .reduce(new BBox(new Vec3(MAX, MAX, MAX), new Vec3(-MAX, -MAX, -MAX)), BBox::union);
  }

  private Hit raycastNode(Node node, Ray ray) {
    if (node == null || !node.bbox.intersects(ray)) return null;

    if (node.isLeaf()) {
      // Check for intersection with all objects in the leaf node.
      return node.objects.stream()
        .map(object -> object.raycast(ray))
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
